#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "XuerichNaIHit.hh"

G4Allocator<XuerichNaIHit> XuerichNaIHitAllocator;

XuerichNaIHit::XuerichNaIHit() {}

XuerichNaIHit::~XuerichNaIHit()
{
	if(m_pParticleType) delete m_pParticleType;
	//if(m_pParentType) delete m_pParentType;
	//if(m_pCreatorProcess) delete m_pCreatorProcess;
	//if(m_pDepositingProcess) delete m_pDepositingProcess;
}

XuerichNaIHit::XuerichNaIHit(const XuerichNaIHit &hXuerichNaIHit):G4VHit()
{
	//m_iTrackId 			= hXuerichNaIHit.m_iTrackId;
	//m_iParentId 			= hXuerichNaIHit.m_iParentId;
	m_pParticleType 		= hXuerichNaIHit.m_pParticleType;
	//m_pParentType 		= hXuerichNaIHit.m_pParentType ;
	//m_pCreatorProcess 	= hXuerichNaIHit.m_pCreatorProcess ;
	//m_pDepositingProcess 	= hXuerichNaIHit.m_pDepositingProcess ;
	m_hPosition		 		= hXuerichNaIHit.m_hPosition;
	m_dEnergyDeposited		= hXuerichNaIHit.m_dEnergyDeposited;
	m_dKineticEnergy 		= hXuerichNaIHit.m_dKineticEnergy ;
	m_dTime		 			= hXuerichNaIHit.m_dTime;
}

const XuerichNaIHit &
XuerichNaIHit::operator=(const XuerichNaIHit &hXuerichNaIHit)
{
	//m_iTrackId 			= hXuerichNaIHit.m_iTrackId;
	//m_iParentId 			= hXuerichNaIHit.m_iParentId;
	m_pParticleType 		= hXuerichNaIHit.m_pParticleType;
	//m_pParentType 		= hXuerichNaIHit.m_pParentType ;
	//m_pCreatorProcess 	= hXuerichNaIHit.m_pCreatorProcess ;
	//m_pDepositingProcess 	= hXuerichNaIHit.m_pDepositingProcess ;
	m_hPosition 			= hXuerichNaIHit.m_hPosition;
	m_dEnergyDeposited 		= hXuerichNaIHit.m_dEnergyDeposited;
	m_dKineticEnergy 		= hXuerichNaIHit.m_dKineticEnergy ;
	m_dTime 				= hXuerichNaIHit.m_dTime;
	
	return *this;
}

G4int
XuerichNaIHit::operator==(const XuerichNaIHit &hXuerichNaIHit) const
{
	return ((this == &hXuerichNaIHit) ? (1) : (0));
}

void XuerichNaIHit::Draw()
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

void XuerichNaIHit::Print()
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

