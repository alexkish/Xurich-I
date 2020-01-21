#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "XuerichLeadHit.hh"

G4Allocator<XuerichLeadHit> XuerichLeadHitAllocator;

XuerichLeadHit::XuerichLeadHit() {}

XuerichLeadHit::~XuerichLeadHit()
{
	if(m_pParticleType) delete m_pParticleType;
	//if(m_pParentType) delete m_pParentType;
	//if(m_pCreatorProcess) delete m_pCreatorProcess;
	//if(m_pDepositingProcess) delete m_pDepositingProcess;
}

XuerichLeadHit::XuerichLeadHit(const XuerichLeadHit &hXuerichLeadHit):G4VHit()
{
	//m_iTrackId = hXuerichLeadHit.m_iTrackId;
	//m_iParentId = hXuerichLeadHit.m_iParentId;
	m_pParticleType = hXuerichLeadHit.m_pParticleType;
	//m_pParentType = hXuerichLeadHit.m_pParentType ;
	//m_pCreatorProcess = hXuerichLeadHit.m_pCreatorProcess ;
	//m_pDepositingProcess = hXuerichLeadHit.m_pDepositingProcess ;
	m_hPosition = hXuerichLeadHit.m_hPosition;
	m_dEnergyDeposited = hXuerichLeadHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXuerichLeadHit.m_dKineticEnergy ;
	m_dTime = hXuerichLeadHit.m_dTime;
}

const XuerichLeadHit &
XuerichLeadHit::operator=(const XuerichLeadHit &hXuerichLeadHit)
{
	//m_iTrackId = hXuerichLeadHit.m_iTrackId;
	//m_iParentId = hXuerichLeadHit.m_iParentId;
	m_pParticleType = hXuerichLeadHit.m_pParticleType;
	//m_pParentType = hXuerichLeadHit.m_pParentType ;
	//m_pCreatorProcess = hXuerichLeadHit.m_pCreatorProcess ;
	//m_pDepositingProcess = hXuerichLeadHit.m_pDepositingProcess ;
	m_hPosition = hXuerichLeadHit.m_hPosition;
	m_dEnergyDeposited = hXuerichLeadHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXuerichLeadHit.m_dKineticEnergy ;
	m_dTime = hXuerichLeadHit.m_dTime;
	
	return *this;
}

G4int
XuerichLeadHit::operator==(const XuerichLeadHit &hXuerichLeadHit) const
{
	return ((this == &hXuerichLeadHit) ? (1) : (0));
}

void XuerichLeadHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
	
	if(pVVisManager)
	{
		G4Circle hCircle(m_hPosition);
		G4Colour hColour(1.000, 0.973, 0.184);
		G4VisAttributes hVisAttributes(hColour);
		
		hCircle.SetScreenSize(0.1);
		hCircle.SetFillStyle(G4Circle::filled);
		hCircle.SetVisAttributes(hVisAttributes);
		pVVisManager->Draw(hCircle);
	}
}

void XuerichLeadHit::Print()
{
	G4cout << "-------------------- LXe hit --------------------" 
		//<< "Id: " << m_iTrackId
		//<< " Particle: " << *m_pParticleType
		//<< " ParentId: " << m_iParentId
		//<< " ParentType: " << *m_pParentType << G4endl
		//<< "CreatorProcess: " << *m_pCreatorProcess
		//<< " DepositingProcess: " << *m_pDepositingProcess << G4endl
		<< "Position: " << m_hPosition.x()/mm
		<< " " << m_hPosition.y()/mm
		<< " " << m_hPosition.z()/mm
		<< " mm" << G4endl
		<< "EnergyDeposited: " << m_dEnergyDeposited/keV << " keV"
		<< " KineticEnergyLeft: " << m_dKineticEnergy/keV << " keV"
		<< " Time: " << m_dTime/s << " s" << G4endl;
}

