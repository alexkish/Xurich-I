#include <G4SDManager.hh>
#include <G4Run.hh>
#include <G4Event.hh>
#include <G4HCofThisEvent.hh>
#include <G4Step.hh>

#include <numeric>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>

#include "XuerichDetectorConstruction.hh"
#include "XuerichPrimaryGeneratorAction.hh"
#include "XuerichEventData.hh"

#include "XuerichLXeHit.hh"
#include "XuerichPmtHit.hh"
#include "XuerichNaIHit.hh"
#include "XuerichTeflonHit.hh"
#include "XuerichDeadXeHit.hh"
#include "XuerichGXeHit.hh"
#include "XuerichSteelHit.hh"
#include "XuerichAlHit.hh"
#include "XuerichLeadHit.hh"
#include "XuerichRestHit.hh"


#include "XuerichAnalysisManager.hh"

XuerichAnalysisManager::XuerichAnalysisManager(XuerichPrimaryGeneratorAction *pPrimaryGeneratorAction)
{
	m_iLXeHitsCollectionID 		= -1;
	m_iPmtHitsCollectionID 		= -1;
	m_iNaIHitsCollectionID 		= -1;
	m_iTeflonHitsCollectionID 	= -1;
	m_iDeadXeHitsCollectionID 	= -1;
	m_iGXeHitsCollectionID 		= -1;
	m_iSteelHitsCollectionID 	= -1;
	m_iAlHitsCollectionID 		= -1;
	m_iLeadHitsCollectionID 	= -1;
	m_iRestHitsCollectionID 	= -1;

	m_hDataFilename = "events.root";

	m_pPrimaryGeneratorAction = pPrimaryGeneratorAction;

	m_pEventData = new XuerichEventData();
}

XuerichAnalysisManager::~XuerichAnalysisManager()
{
}

void
XuerichAnalysisManager::BeginOfRun(const G4Run *pRun)
{
	m_pTreeFile = new TFile(m_hDataFilename.c_str(), "RECREATE", "File containing event data for Xuerich");
	m_pTree = new TTree("t1", "Tree containing event data for Xuerich");

	gROOT->ProcessLine("#include <vector>");

	m_pTree->Branch("eventid", &m_pEventData->m_iEventId, "eventid/I");
	m_pTree->Branch("ntpmthits", &m_pEventData->m_iNbTopPmtHits, "ntpmthits/I");
	m_pTree->Branch("nbpmthits", &m_pEventData->m_iNbBottomPmtHits, "nbpmthits/I");
	m_pTree->Branch("pmthits", "vector<int>", &m_pEventData->m_pPmtHits);
	m_pTree->Branch("nsteps", &m_pEventData->m_iNbSteps, "nsteps/I");

	//m_pTree->Branch("etotER", &m_pEventData->m_fTotalEnergyDepositedER, "etotER/F");
	//m_pTree->Branch("etotNR", &m_pEventData->m_fTotalEnergyDepositedNR, "etotNR/F");
	m_pTree->Branch("etot", &m_pEventData->m_fTotalEnergyDeposited, "etot/F");
	
	m_pTree->Branch("trackid", "vector<int>", &m_pEventData->m_pTrackId);
	m_pTree->Branch("type", "vector<string>", &m_pEventData->m_pParticleType);
	m_pTree->Branch("parentid", "vector<int>", &m_pEventData->m_pParentId);
	m_pTree->Branch("parenttype", "vector<string>", &m_pEventData->m_pParentType);
	m_pTree->Branch("creaproc", "vector<string>", &m_pEventData->m_pCreatorProcess);
	m_pTree->Branch("edproc", "vector<string>", &m_pEventData->m_pDepositingProcess);
	m_pTree->Branch("xp", "vector<float>", &m_pEventData->m_pX);
	m_pTree->Branch("yp", "vector<float>", &m_pEventData->m_pY);
	m_pTree->Branch("zp", "vector<float>", &m_pEventData->m_pZ);
	m_pTree->Branch("ed", "vector<float>", &m_pEventData->m_pEnergyDeposited);
	m_pTree->Branch("time", "vector<float>", &m_pEventData->m_pTime);

	m_pTree->Branch("pre_kinE", "vector<double>", &m_pEventData->m_pPreKinNrg);
	m_pTree->Branch("post_kinE", "vector<double>", &m_pEventData->m_pPostKinNrg);

	m_pTree->Branch("type_pri", "vector<string>", &m_pEventData->m_pPrimaryParticleType);
	m_pTree->Branch("xp_pri", &m_pEventData->m_fPrimaryX, "xp_pri/F");
	m_pTree->Branch("yp_pri", &m_pEventData->m_fPrimaryY, "yp_pri/F");
	m_pTree->Branch("zp_pri", &m_pEventData->m_fPrimaryZ, "zp_pri/F");

	m_pTree->Branch("NaI_nsteps", 	&m_pEventData->m_iNbStepsNaI, "NaI_nsteps/I");
	m_pTree->Branch("NaI_ed", 		"vector<float>", &m_pEventData->m_pEnergyDepositedNaI);
	//m_pTree->Branch("NaI_etotER", 	&m_pEventData->m_fTotalEnergyDepositedNaIER, "NaI_etotER/F");
	//m_pTree->Branch("NaI_etotNR", 	&m_pEventData->m_fTotalEnergyDepositedNaINR, "NaI_etotNR/F");
	m_pTree->Branch("NaI_etot", 	&m_pEventData->m_fTotalEnergyDepositedNaI, "NaI_etot/F");
	m_pTree->Branch("NaI_time", 	"vector<float>", &m_pEventData->m_pTimeNaI);
	m_pTree->Branch("NaI_xp", 		"vector<float>", &m_pEventData->m_pXNaI);
	m_pTree->Branch("NaI_yp", 		"vector<float>", &m_pEventData->m_pYNaI);
	m_pTree->Branch("NaI_zp", 		"vector<float>", &m_pEventData->m_pZNaI);

	m_pTree->Branch("Teflon_nsteps", &m_pEventData->m_iNbStepsTeflon, "Teflon_nsteps/I");
	m_pTree->Branch("Teflon_ed", "vector<float>", &m_pEventData->m_pEnergyDepositedTeflon);
	m_pTree->Branch("Teflon_etot", &m_pEventData->m_fTotalEnergyDepositedTeflon, "Teflon_etot/F");
	m_pTree->Branch("Teflon_time", "vector<float>", &m_pEventData->m_pTimeTeflon);

	m_pTree->Branch("DeadXe_nsteps", &m_pEventData->m_iNbStepsDeadXe, "DeadXe_nsteps/I");
	m_pTree->Branch("DeadXe_ed", "vector<float>", &m_pEventData->m_pEnergyDepositedDeadXe);
	m_pTree->Branch("DeadXe_etot", &m_pEventData->m_fTotalEnergyDepositedDeadXe, "DeadXe_etot/F");
	m_pTree->Branch("DeadXe_time", "vector<float>", &m_pEventData->m_pTimeDeadXe);

	m_pTree->Branch("GXe_nsteps", &m_pEventData->m_iNbStepsGXe, "GXe_nsteps/I");
	m_pTree->Branch("GXe_ed", "vector<float>", &m_pEventData->m_pEnergyDepositedGXe);
	m_pTree->Branch("GXe_etot", &m_pEventData->m_fTotalEnergyDepositedGXe, "GXe_etot/F");
	m_pTree->Branch("GXe_time", "vector<float>", &m_pEventData->m_pTimeGXe);

	m_pTree->Branch("Steel_nsteps", &m_pEventData->m_iNbStepsSteel, "Steel_nsteps/I");
	m_pTree->Branch("Steel_ed", "vector<float>", &m_pEventData->m_pEnergyDepositedSteel);
	m_pTree->Branch("Steel_etot", &m_pEventData->m_fTotalEnergyDepositedSteel, "Steel_etot/F");
	m_pTree->Branch("Steel_time", "vector<float>", &m_pEventData->m_pTimeSteel);

	m_pTree->Branch("Al_nsteps", &m_pEventData->m_iNbStepsAl, "Al_nsteps/I");
	m_pTree->Branch("Al_ed", "vector<float>", &m_pEventData->m_pEnergyDepositedAl);
	m_pTree->Branch("Al_etot", &m_pEventData->m_fTotalEnergyDepositedAl, "Al_etot/F");
	m_pTree->Branch("Al_time", "vector<float>", &m_pEventData->m_pTimeAl);

	m_pTree->Branch("Lead_nsteps", &m_pEventData->m_iNbStepsLead, "Lead_nsteps/I");
	m_pTree->Branch("Lead_ed", "vector<float>", &m_pEventData->m_pEnergyDepositedLead);
	m_pTree->Branch("Lead_etot", &m_pEventData->m_fTotalEnergyDepositedLead, "Lead_etot/F");
	m_pTree->Branch("Lead_time", "vector<float>", &m_pEventData->m_pTimeLead);

	m_pTree->Branch("Rest_nsteps", &m_pEventData->m_iNbStepsRest, "Rest_nsteps/I");
	m_pTree->Branch("Rest_ed", "vector<float>", &m_pEventData->m_pEnergyDepositedRest);
	m_pTree->Branch("Rest_etot", &m_pEventData->m_fTotalEnergyDepositedRest, "Rest_etot/F");
	m_pTree->Branch("Rest_time", "vector<float>", &m_pEventData->m_pTimeRest);

	m_pTree->SetMaxTreeSize(10e9);
	m_pTree->AutoSave();

	m_pNbEventsToSimulateParameter = new TParameter<int>("nbevents", m_iNbEventsToSimulate);
	m_pNbEventsToSimulateParameter->Write();
}

void
XuerichAnalysisManager::EndOfRun(const G4Run *pRun)
{
	m_pTreeFile->Write();
	m_pTreeFile->Close();
}

void
XuerichAnalysisManager::BeginOfEvent(const G4Event *pEvent)
{
	if(m_iLXeHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iLXeHitsCollectionID 		= pSDManager->GetCollectionID("LXeHitsCollection");
	} 

	if(m_iPmtHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iPmtHitsCollectionID 		= pSDManager->GetCollectionID("PmtHitsCollection");
	}

	if(m_iNaIHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iNaIHitsCollectionID 		= pSDManager->GetCollectionID("NaIHitsCollection");
	} 

	if(m_iTeflonHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iTeflonHitsCollectionID 	= pSDManager->GetCollectionID("TeflonHitsCollection");
	} 

	if(m_iDeadXeHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iDeadXeHitsCollectionID 	= pSDManager->GetCollectionID("DeadXeHitsCollection");
	} 

	if(m_iGXeHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iGXeHitsCollectionID 		= pSDManager->GetCollectionID("GXeHitsCollection");
	} 

	if(m_iSteelHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iSteelHitsCollectionID 	= pSDManager->GetCollectionID("SteelHitsCollection");
	} 

	if(m_iAlHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iAlHitsCollectionID 		= pSDManager->GetCollectionID("AlHitsCollection");
	} 

	if(m_iLeadHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iLeadHitsCollectionID 	= pSDManager->GetCollectionID("LeadHitsCollection");
	} 

	if(m_iRestHitsCollectionID == -1)
	{
		G4SDManager *pSDManager 	= G4SDManager::GetSDMpointer();
		m_iRestHitsCollectionID 	= pSDManager->GetCollectionID("RestHitsCollection");
	} 

}

void
XuerichAnalysisManager::EndOfEvent(const G4Event *pEvent)
{
	G4HCofThisEvent* pHCofThisEvent = pEvent->GetHCofThisEvent();
	XuerichLXeHitsCollection	* pLXeHitsCollection 	= 0;
	XuerichPmtHitsCollection	* pPmtHitsCollection 	= 0;
	XuerichNaIHitsCollection	* pNaIHitsCollection 	= 0;
	XuerichTeflonHitsCollection	* pTeflonHitsCollection = 0;
	XuerichDeadXeHitsCollection	* pDeadXeHitsCollection = 0;
	XuerichGXeHitsCollection	* pGXeHitsCollection 	= 0;
	XuerichSteelHitsCollection	* pSteelHitsCollection 	= 0;
	XuerichAlHitsCollection		* pAlHitsCollection 	= 0;
	XuerichLeadHitsCollection	* pLeadHitsCollection 	= 0;
	XuerichRestHitsCollection	* pRestHitsCollection 	= 0;

	G4int iNbLXeHits 	= 0;
	G4int iNbPmtHits 	= 0;
	G4int iNbNaIHits 	= 0;
	G4int iNbTeflonHits = 0;
	G4int iNbDeadXeHits = 0;
	G4int iNbGXeHits 	= 0;
	G4int iNbSteelHits 	= 0;
	G4int iNbAlHits 	= 0;
	G4int iNbLeadHits 	= 0;
	G4int iNbRestHits 	= 0;
	
	if(pHCofThisEvent)
	{
		if(m_iLXeHitsCollectionID != -1)
		{
			pLXeHitsCollection = (XuerichLXeHitsCollection *)(pHCofThisEvent->GetHC(m_iLXeHitsCollectionID));
			iNbLXeHits = (pLXeHitsCollection)?(pLXeHitsCollection->entries()):(0);
		}

		if(m_iPmtHitsCollectionID != -1)
		{
			pPmtHitsCollection = (XuerichPmtHitsCollection *)(pHCofThisEvent->GetHC(m_iPmtHitsCollectionID));
			iNbPmtHits = (pPmtHitsCollection)?(pPmtHitsCollection->entries()):(0);
		}

		if(m_iNaIHitsCollectionID != -1)
		{
			pNaIHitsCollection = (XuerichNaIHitsCollection *)(pHCofThisEvent->GetHC(m_iNaIHitsCollectionID));
			iNbNaIHits = (pNaIHitsCollection)?(pNaIHitsCollection->entries()):(0);
		}

		if(m_iTeflonHitsCollectionID != -1)
		{
			pTeflonHitsCollection = (XuerichTeflonHitsCollection *)(pHCofThisEvent->GetHC(m_iTeflonHitsCollectionID));
			iNbTeflonHits = (pTeflonHitsCollection)?(pTeflonHitsCollection->entries()):(0);
		}

		if(m_iDeadXeHitsCollectionID != -1)
		{
			pDeadXeHitsCollection = (XuerichDeadXeHitsCollection *)(pHCofThisEvent->GetHC(m_iDeadXeHitsCollectionID));
			iNbDeadXeHits = (pDeadXeHitsCollection)?(pDeadXeHitsCollection->entries()):(0);
		}

		if(m_iGXeHitsCollectionID != -1)
		{
			pGXeHitsCollection = (XuerichGXeHitsCollection *)(pHCofThisEvent->GetHC(m_iGXeHitsCollectionID));
			iNbGXeHits = (pGXeHitsCollection)?(pGXeHitsCollection->entries()):(0);
		}

		if(m_iSteelHitsCollectionID != -1)
		{
			pSteelHitsCollection = (XuerichSteelHitsCollection *)(pHCofThisEvent->GetHC(m_iSteelHitsCollectionID));
			iNbSteelHits = (pSteelHitsCollection)?(pSteelHitsCollection->entries()):(0);
		}

		if(m_iAlHitsCollectionID != -1)
		{
			pAlHitsCollection = (XuerichAlHitsCollection *)(pHCofThisEvent->GetHC(m_iAlHitsCollectionID));
			iNbAlHits = (pAlHitsCollection)?(pAlHitsCollection->entries()):(0);
		}

		if(m_iLeadHitsCollectionID != -1)
		{
			pLeadHitsCollection = (XuerichLeadHitsCollection *)(pHCofThisEvent->GetHC(m_iLeadHitsCollectionID));
			iNbLeadHits = (pLeadHitsCollection)?(pLeadHitsCollection->entries()):(0);
		}

		if(m_iRestHitsCollectionID != -1)
		{
			pRestHitsCollection = (XuerichRestHitsCollection *)(pHCofThisEvent->GetHC(m_iRestHitsCollectionID));
			iNbRestHits = (pRestHitsCollection)?(pRestHitsCollection->entries()):(0);
		}

	}

	//if(iNbLXeHits)
	if(iNbLXeHits>0 && iNbNaIHits>0)
	{
		m_pEventData->m_iEventId = pEvent->GetEventID();

		m_pEventData->m_pPrimaryParticleType->push_back(m_pPrimaryGeneratorAction->GetParticleTypeOfPrimary());

		m_pEventData->m_fPrimaryX = m_pPrimaryGeneratorAction->GetPositionOfPrimary().x();
		m_pEventData->m_fPrimaryY = m_pPrimaryGeneratorAction->GetPositionOfPrimary().y();
		m_pEventData->m_fPrimaryZ = m_pPrimaryGeneratorAction->GetPositionOfPrimary().z();

		G4int iNbSteps 			= 0;
		G4int iNbStepsNaI 		= 0;
		G4int iNbStepsTeflon 	= 0;
		G4int iNbStepsDeadXe 	= 0;
		G4int iNbStepsGXe 		= 0;
		G4int iNbStepsSteel 	= 0;
		G4int iNbStepsAl 		= 0;
		G4int iNbStepsLead 		= 0;
		G4int iNbStepsRest 		= 0;

		//G4float fTotalEnergyDepositedER 	= 0.;
		//G4float fTotalEnergyDepositedNR 	= 0.;
		//G4float fTotalEnergyDepositedNaIER 	= 0.;
		//G4float fTotalEnergyDepositedNaINR 	= 0.;
		G4float fTotalEnergyDeposited 		= 0.;
		G4float fTotalEnergyDepositedNaI 	= 0.;
		G4float fTotalEnergyDepositedTeflon = 0.;
		G4float fTotalEnergyDepositedDeadXe = 0.;
		G4float fTotalEnergyDepositedGXe 	= 0.;
		G4float fTotalEnergyDepositedSteel 	= 0.;
		G4float fTotalEnergyDepositedAl 	= 0.;
		G4float fTotalEnergyDepositedLead 	= 0.;
		G4float fTotalEnergyDepositedRest 	= 0.;

		// LXe hits
		for(G4int i=0; i<iNbLXeHits; i++)
		{
			XuerichLXeHit *pHit = (*pLXeHitsCollection)[i];

		/*	if(pHit->GetParticleType()=="gamma" || pHit->GetParticleType()=="e-" || pHit->GetParticleType()=="e+")
			{
				m_pEventData->m_pTrackId->push_back(pHit->GetTrackId());
				m_pEventData->m_pParentId->push_back(pHit->GetParentId());

				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());
				m_pEventData->m_pParentType->push_back(pHit->GetParentType());
				m_pEventData->m_pCreatorProcess->push_back(pHit->GetCreatorProcess());
				m_pEventData->m_pDepositingProcess->push_back(pHit->GetDepositingProcess());

				m_pEventData->m_pX->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pY->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZ->push_back(pHit->GetPosition().z()/mm);

				fTotalEnergyDepositedER += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDeposited->push_back(pHit->GetEnergyDeposited()/keV);

				m_pEventData->m_pKineticEnergy->push_back(pHit->GetKineticEnergy()/keV);
				m_pEventData->m_pTime->push_back(pHit->GetTime()/second);

				iNbSteps++;
			}

			if(pHit->GetParticleType()!="gamma" && pHit->GetParticleType()!="e-" && pHit->GetParticleType()!="e+")
			{
				m_pEventData->m_pTrackId->push_back(pHit->GetTrackId());
				m_pEventData->m_pParentId->push_back(pHit->GetParentId());

				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());
				m_pEventData->m_pParentType->push_back(pHit->GetParentType());
				m_pEventData->m_pCreatorProcess->push_back(pHit->GetCreatorProcess());
				m_pEventData->m_pDepositingProcess->push_back(pHit->GetDepositingProcess());

				m_pEventData->m_pX->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pY->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZ->push_back(pHit->GetPosition().z()/mm);

				fTotalEnergyDepositedNR += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDeposited->push_back(pHit->GetEnergyDeposited()/keV);

				m_pEventData->m_pKineticEnergy->push_back(pHit->GetKineticEnergy()/keV);
				m_pEventData->m_pTime->push_back(pHit->GetTime()/second);

				iNbSteps++;
			}
		*/
			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pTrackId->push_back(pHit->GetTrackId());
				m_pEventData->m_pParentId->push_back(pHit->GetParentId());

				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());
				m_pEventData->m_pParentType->push_back(pHit->GetParentType());
				m_pEventData->m_pCreatorProcess->push_back(pHit->GetCreatorProcess());
				m_pEventData->m_pDepositingProcess->push_back(pHit->GetDepositingProcess());

				m_pEventData->m_pX->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pY->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZ->push_back(pHit->GetPosition().z()/mm);

				m_pEventData->m_pPreKinNrg->push_back(pHit->GetPreKinEnergy()/keV);
				m_pEventData->m_pPostKinNrg->push_back(pHit->GetPostKinEnergy()/keV);

				fTotalEnergyDeposited += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDeposited->push_back(pHit->GetEnergyDeposited()/keV);

				m_pEventData->m_pKineticEnergy->push_back(pHit->GetKineticEnergy()/keV);
				
				m_pEventData->m_pTime->push_back(pHit->GetTime()/second);

				iNbSteps++;
			}
		};

		// NaI scintillator hits
		for(G4int i=0; i<iNbNaIHits; i++)
		{
			XuerichNaIHit *pHit = (*pNaIHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			//if(pHit->GetParticleType()=="gamma" || pHit->GetParticleType()=="e-" || pHit->GetParticleType()=="e+")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedNaI += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedNaI->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeNaI->push_back(pHit->GetTime()/second);

				m_pEventData->m_pXNaI->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pYNaI->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZNaI->push_back(pHit->GetPosition().z()/mm);

				iNbStepsNaI++;
			}

		/*	if(pHit->GetParticleType()!="gamma" && pHit->GetParticleType()!="e-" && pHit->GetParticleType()!="e+")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedNaINR += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedNaI->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeNaI->push_back(pHit->GetTime()/second);

				m_pEventData->m_pXNaI->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pYNaI->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZNaI->push_back(pHit->GetPosition().z()/mm);

				iNbStepsNaI++;
			}
		*/
		};

		// teflon hits
		for(G4int i=0; i<iNbTeflonHits; i++)
		{
			XuerichTeflonHit *pHit = (*pTeflonHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedTeflon += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedTeflon->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeTeflon->push_back(pHit->GetTime()/second);

				iNbStepsTeflon++;
			}
		};

		// dead Xe hits
		for(G4int i=0; i<iNbDeadXeHits; i++)
		{
			XuerichDeadXeHit *pHit = (*pDeadXeHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedDeadXe += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedDeadXe->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeDeadXe->push_back(pHit->GetTime()/second);

				iNbStepsDeadXe++;
			}
		};

		// gas GXe hits
		for(G4int i=0; i<iNbGXeHits; i++)
		{
			XuerichGXeHit *pHit = (*pGXeHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedGXe += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedGXe->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeGXe->push_back(pHit->GetTime()/second);

				iNbStepsGXe++;
			}
		};

		// steel hits
		for(G4int i=0; i<iNbSteelHits; i++)
		{
			XuerichSteelHit *pHit = (*pSteelHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedSteel += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedSteel->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeSteel->push_back(pHit->GetTime()/second);

				iNbStepsSteel++;
			}
		};

		// Aluminum hits
		for(G4int i=0; i<iNbAlHits; i++)
		{
			XuerichAlHit *pHit = (*pAlHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedAl += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedAl->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeAl->push_back(pHit->GetTime()/second);

				iNbStepsAl++;
			}
		};

		// Lead hits
		for(G4int i=0; i<iNbLeadHits; i++)
		{
			XuerichLeadHit *pHit = (*pLeadHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedLead += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedLead->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeLead->push_back(pHit->GetTime()/second);

				iNbStepsLead++;
			}
		};

		// hits in the rest of materials
		for(G4int i=0; i<iNbRestHits; i++)
		{
			XuerichRestHit *pHit = (*pRestHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());

				fTotalEnergyDepositedRest += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDepositedRest->push_back(pHit->GetEnergyDeposited()/keV);
				m_pEventData->m_pTimeRest->push_back(pHit->GetTime()/second);

				iNbStepsRest++;
			}
		};

		m_pEventData->m_iNbSteps 		= iNbSteps;
		m_pEventData->m_iNbStepsNaI 	= iNbStepsNaI;
		m_pEventData->m_iNbStepsTeflon 	= iNbStepsTeflon;
		m_pEventData->m_iNbStepsDeadXe 	= iNbStepsDeadXe;
		m_pEventData->m_iNbStepsGXe 	= iNbStepsGXe;
		m_pEventData->m_iNbStepsSteel 	= iNbStepsSteel;
		m_pEventData->m_iNbStepsAl 		= iNbStepsAl;
		m_pEventData->m_iNbStepsLead	= iNbStepsLead;
		m_pEventData->m_iNbStepsRest 	= iNbStepsRest;

		//m_pEventData->m_fTotalEnergyDepositedER 		= fTotalEnergyDepositedER;
		//m_pEventData->m_fTotalEnergyDepositedNR 		= fTotalEnergyDepositedNR;
		//m_pEventData->m_fTotalEnergyDepositedNaIER  	= fTotalEnergyDepositedNaIER;
		//m_pEventData->m_fTotalEnergyDepositedNaINR  	= fTotalEnergyDepositedNaINR;
		m_pEventData->m_fTotalEnergyDeposited 			= fTotalEnergyDeposited;
		m_pEventData->m_fTotalEnergyDepositedNaI 		= fTotalEnergyDepositedNaI;
		m_pEventData->m_fTotalEnergyDepositedTeflon 	= fTotalEnergyDepositedTeflon;
		m_pEventData->m_fTotalEnergyDepositedDeadXe 	= fTotalEnergyDepositedDeadXe;
		m_pEventData->m_fTotalEnergyDepositedGXe 		= fTotalEnergyDepositedGXe;
		m_pEventData->m_fTotalEnergyDepositedSteel 		= fTotalEnergyDepositedSteel;
		m_pEventData->m_fTotalEnergyDepositedAl 		= fTotalEnergyDepositedAl;
		m_pEventData->m_fTotalEnergyDepositedLead 		= fTotalEnergyDepositedLead;
		m_pEventData->m_fTotalEnergyDepositedRest 		= fTotalEnergyDepositedRest;

/*		const G4int NbPmts = 2;
		//m_pEventData->m_pPmtHits->resize(iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts+iNbBottomVetoPmts, 0);
		// Pmt hits
		for(G4int i=0; i<iNbPmtHits; i++)
			(*(m_pEventData->m_pPmtHits))[(*pPmtHitsCollection)[i]->GetPmtNb()]++;

		m_pEventData->m_iNbTopPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin(), m_pEventData->m_pPmtHits->begin()+iNbTopPmts, 0);
		m_pEventData->m_iNbBottomPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts, m_pEventData->m_pPmtHits->end(), 0);
*/
		if(fTotalEnergyDeposited > 0.)
		//if(fTotalEnergyDepositedNaI > 0.)
		//if(fTotalEnergyDeposited > 0. || fTotalEnergyDepositedNaI > 0.)
		m_pTree->Fill();


		m_pEventData->Clear();
	}
}

void
XuerichAnalysisManager::Step(const G4Step *pStep)
{
}
	
