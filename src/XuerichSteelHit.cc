#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "XuerichSteelHit.hh"

G4Allocator<XuerichSteelHit> XuerichSteelHitAllocator;

XuerichSteelHit::XuerichSteelHit() {}

XuerichSteelHit::~XuerichSteelHit()
{
	if(m_pParticleType) delete m_pParticleType;
	//if(m_pParentType) delete m_pParentType;
	//if(m_pCreatorProcess) delete m_pCreatorProcess;
	//if(m_pDepositingProcess) delete m_pDepositingProcess;
}

XuerichSteelHit::XuerichSteelHit(const XuerichSteelHit &hXuerichSteelHit):G4VHit()
{
	//m_iTrackId = hXuerichSteelHit.m_iTrackId;
	//m_iParentId = hXuerichSteelHit.m_iParentId;
	m_pParticleType = hXuerichSteelHit.m_pParticleType;
	//m_pParentType = hXuerichSteelHit.m_pParentType ;
	//m_pCreatorProcess = hXuerichSteelHit.m_pCreatorProcess ;
	//m_pDepositingProcess = hXuerichSteelHit.m_pDepositingProcess ;
	m_hPosition = hXuerichSteelHit.m_hPosition;
	m_dEnergyDeposited = hXuerichSteelHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXuerichSteelHit.m_dKineticEnergy ;
	m_dTime = hXuerichSteelHit.m_dTime;
}

const XuerichSteelHit &
XuerichSteelHit::operator=(const XuerichSteelHit &hXuerichSteelHit)
{
	//m_iTrackId = hXuerichSteelHit.m_iTrackId;
	//m_iParentId = hXuerichSteelHit.m_iParentId;
	m_pParticleType = hXuerichSteelHit.m_pParticleType;
	//m_pParentType = hXuerichSteelHit.m_pParentType ;
	//m_pCreatorProcess = hXuerichSteelHit.m_pCreatorProcess ;
	//m_pDepositingProcess = hXuerichSteelHit.m_pDepositingProcess ;
	m_hPosition = hXuerichSteelHit.m_hPosition;
	m_dEnergyDeposited = hXuerichSteelHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXuerichSteelHit.m_dKineticEnergy ;
	m_dTime = hXuerichSteelHit.m_dTime;
	
	return *this;
}

G4int
XuerichSteelHit::operator==(const XuerichSteelHit &hXuerichSteelHit) const
{
	return ((this == &hXuerichSteelHit) ? (1) : (0));
}

void XuerichSteelHit::Draw()
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

void XuerichSteelHit::Print()
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

