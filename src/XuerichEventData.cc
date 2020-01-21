#include "XuerichEventData.hh"

XuerichEventData::XuerichEventData()
{
	m_iEventId 			= 0;
	m_iNbStepsNaI 		= 0;
	m_iNbStepsTeflon 	= 0;
	m_iNbStepsDeadXe 	= 0;
	m_iNbStepsGXe 		= 0;
	m_iNbStepsSteel 	= 0;
	m_iNbStepsAl 		= 0;
	m_iNbStepsLead 		= 0;
	m_iNbStepsRest 		= 0;

	m_iNbTopPmtHits 	= 0;
	m_iNbBottomPmtHits 	= 0;
	m_pPmtHits 			= new vector<int>;

	//m_fTotalEnergyDepositedER 	= 0.;
	//m_fTotalEnergyDepositedNR 	= 0.;
	//m_fTotalEnergyDepositedNaIER 	= 0.;
	//m_fTotalEnergyDepositedNaINR 	= 0.;

	m_fTotalEnergyDeposited 		= 0.;
	m_fTotalEnergyDepositedNaI 		= 0.;
	m_fTotalEnergyDepositedTeflon 	= 0.;
	m_fTotalEnergyDepositedDeadXe 	= 0.;
	m_fTotalEnergyDepositedGXe 		= 0.;
	m_fTotalEnergyDepositedSteel 	= 0.;
	m_fTotalEnergyDepositedAl 		= 0.;
	m_fTotalEnergyDepositedLead 	= 0.;
	m_fTotalEnergyDepositedRest 	= 0.;


	m_pTrackId = new vector<int>;
	m_pParentId = new vector<int>;

	m_pParticleType = new vector<string>;
	m_pParticleTypeNaI = new vector<string>;

	m_pParentType = new vector<string>;
	m_pCreatorProcess = new vector<string>;
	m_pDepositingProcess = new vector<string>;

	m_pX = new vector<float>;
	m_pY = new vector<float>;
	m_pZ = new vector<float>;

	m_pXNaI = new vector<float>;
	m_pYNaI = new vector<float>;
	m_pZNaI = new vector<float>;

	m_pEnergyDeposited 			= new vector<float>;
	m_pEnergyDepositedNaI 		= new vector<float>;
	m_pEnergyDepositedTeflon 	= new vector<float>;
	m_pEnergyDepositedDeadXe 	= new vector<float>;
	m_pEnergyDepositedGXe 		= new vector<float>;
	m_pEnergyDepositedSteel 	= new vector<float>;
	m_pEnergyDepositedAl 		= new vector<float>;
	m_pEnergyDepositedLead 		= new vector<float>;
	m_pEnergyDepositedRest 		= new vector<float>;

	m_pKineticEnergy = new vector<float>;
    m_pPreKinNrg = new vector<double>;
    m_pPostKinNrg = new vector<double>;

	m_pTime 		= new vector<float>;
	m_pTimeNaI 		= new vector<float>;
	m_pTimeTeflon 	= new vector<float>;
	m_pTimeDeadXe 	= new vector<float>;
	m_pTimeGXe 		= new vector<float>;
	m_pTimeSteel 	= new vector<float>;
	m_pTimeAl 		= new vector<float>;
	m_pTimeLead 	= new vector<float>;
	m_pTimeRest 	= new vector<float>;

	m_pPrimaryParticleType = new vector<string>;
	m_fPrimaryX = 0.;
	m_fPrimaryY = 0.;
	m_fPrimaryZ = 0.;	
}

XuerichEventData::~XuerichEventData()
{
	delete m_pPmtHits;
	delete m_pTrackId;
	delete m_pParentId;
	
	delete m_pParticleType;
	delete m_pParticleTypeNaI;

	delete m_pParentType;
	delete m_pCreatorProcess;
	delete m_pDepositingProcess;
	delete m_pKineticEnergy;
    delete m_pPreKinNrg;
    delete m_pPostKinNrg;

	delete m_pX;
	delete m_pY;
	delete m_pZ;

	delete m_pXNaI;
	delete m_pYNaI;
	delete m_pZNaI;

	delete m_pEnergyDeposited;
	delete m_pEnergyDepositedNaI;
	delete m_pEnergyDepositedTeflon;
	delete m_pEnergyDepositedDeadXe;
	delete m_pEnergyDepositedGXe;
	delete m_pEnergyDepositedSteel;
	delete m_pEnergyDepositedAl;
	delete m_pEnergyDepositedLead;
	delete m_pEnergyDepositedRest;

	delete m_pTime;
	delete m_pTimeNaI;
	delete m_pTimeTeflon;
	delete m_pTimeDeadXe;
	delete m_pTimeGXe;
	delete m_pTimeSteel;
	delete m_pTimeAl;
	delete m_pTimeRest;

	delete m_pPrimaryParticleType;
}

void
XuerichEventData::Clear()
{
	m_iEventId 			= 0;
	m_iNbSteps 			= 0;
	m_iNbStepsNaI 		= 0;
	m_iNbStepsTeflon 	= 0;
	m_iNbStepsDeadXe 	= 0;
	m_iNbStepsGXe 		= 0;
	m_iNbStepsSteel 	= 0;
	m_iNbStepsAl 		= 0;
	m_iNbStepsLead 		= 0;
	m_iNbStepsRest 		= 0;


	m_iNbTopPmtHits = 0;
	m_iNbBottomPmtHits = 0;

	m_pPmtHits->clear();

	//m_fTotalEnergyDepositedER 	= 0.0;
	//m_fTotalEnergyDepositedNR 	= 0.0;
	//m_fTotalEnergyDepositedNaIER 	= 0.0;
	//m_fTotalEnergyDepositedNaINR 	= 0.0;
	m_fTotalEnergyDeposited 		= 0.0;
	m_fTotalEnergyDepositedNaI 		= 0.0;
	m_fTotalEnergyDepositedTeflon 	= 0.0;
	m_fTotalEnergyDepositedDeadXe 	= 0.0;
	m_fTotalEnergyDepositedGXe 		= 0.0;
	m_fTotalEnergyDepositedSteel 	= 0.0;
	m_fTotalEnergyDepositedAl 		= 0.0;
	m_fTotalEnergyDepositedLead 	= 0.0;
	m_fTotalEnergyDepositedRest 	= 0.0;

	m_pTrackId->clear();
	m_pParentId->clear();

	m_pParticleType->clear();
	m_pParticleTypeNaI->clear();

	m_pParentType->clear();
	m_pCreatorProcess->clear();
	m_pDepositingProcess->clear();

	m_pX->clear();
	m_pY->clear();
	m_pZ->clear();

	m_pXNaI->clear();
	m_pYNaI->clear();
	m_pZNaI->clear();

	m_pEnergyDeposited->clear();
	m_pEnergyDepositedNaI->clear();
	m_pEnergyDepositedTeflon->clear();
	m_pEnergyDepositedDeadXe->clear();
	m_pEnergyDepositedGXe->clear();
	m_pEnergyDepositedSteel->clear();
	m_pEnergyDepositedAl->clear();
	m_pEnergyDepositedLead->clear();
	m_pEnergyDepositedRest->clear();
	m_pKineticEnergy->clear();
    m_pPreKinNrg->clear();
    m_pPostKinNrg->clear();

	m_pTime->clear();
	m_pTimeNaI->clear();
	m_pTimeTeflon->clear();
	m_pTimeDeadXe->clear();
	m_pTimeGXe->clear();
	m_pTimeSteel->clear();
	m_pTimeAl->clear();
	m_pTimeLead->clear();
	m_pTimeRest->clear();

	m_pPrimaryParticleType->clear();
	m_fPrimaryX = 0.;
	m_fPrimaryY = 0.;
	m_fPrimaryZ = 0.;	
}

