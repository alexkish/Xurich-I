#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "XuerichAlHit.hh"

G4Allocator<XuerichAlHit> XuerichAlHitAllocator;

XuerichAlHit::XuerichAlHit() {}

XuerichAlHit::~XuerichAlHit()
{
	if(m_pParticleType) delete m_pParticleType;
	//if(m_pParentType) delete m_pParentType;
	//if(m_pCreatorProcess) delete m_pCreatorProcess;
	//if(m_pDepositingProcess) delete m_pDepositingProcess;
}

XuerichAlHit::XuerichAlHit(const XuerichAlHit &hXuerichAlHit):G4VHit()
{
	//m_iTrackId = hXuerichAlHit.m_iTrackId;
	//m_iParentId = hXuerichAlHit.m_iParentId;
	m_pParticleType = hXuerichAlHit.m_pParticleType;
	//m_pParentType = hXuerichAlHit.m_pParentType ;
	//m_pCreatorProcess = hXuerichAlHit.m_pCreatorProcess ;
	//m_pDepositingProcess = hXuerichAlHit.m_pDepositingProcess ;
	m_hPosition = hXuerichAlHit.m_hPosition;
	m_dEnergyDeposited = hXuerichAlHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXuerichAlHit.m_dKineticEnergy ;
	m_dTime = hXuerichAlHit.m_dTime;
}

const XuerichAlHit &
XuerichAlHit::operator=(const XuerichAlHit &hXuerichAlHit)
{
	//m_iTrackId = hXuerichAlHit.m_iTrackId;
	//m_iParentId = hXuerichAlHit.m_iParentId;
	m_pParticleType = hXuerichAlHit.m_pParticleType;
	//m_pParentType = hXuerichAlHit.m_pParentType ;
	//m_pCreatorProcess = hXuerichAlHit.m_pCreatorProcess ;
	//m_pDepositingProcess = hXuerichAlHit.m_pDepositingProcess ;
	m_hPosition = hXuerichAlHit.m_hPosition;
	m_dEnergyDeposited = hXuerichAlHit.m_dEnergyDeposited;
	m_dKineticEnergy = hXuerichAlHit.m_dKineticEnergy ;
	m_dTime = hXuerichAlHit.m_dTime;
	
	return *this;
}

G4int
XuerichAlHit::operator==(const XuerichAlHit &hXuerichAlHit) const
{
	return ((this == &hXuerichAlHit) ? (1) : (0));
}

void XuerichAlHit::Draw()
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

void XuerichAlHit::Print()
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

