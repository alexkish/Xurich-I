#include <G4PrimaryParticle.hh>
#include <G4Event.hh>
#include <G4TransportationManager.hh>
#include <G4VPhysicalVolume.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4IonTable.hh>
#include <G4Ions.hh>
#include <G4TrackingManager.hh>
#include <G4Track.hh>
#include <Randomize.hh>

#include <sstream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

using std::stringstream;
using std::vector;
using std::ifstream;

#include "XuerichParticleSource.hh"

XuerichParticleSource::XuerichParticleSource()
{
	particleTable = G4ParticleTable::GetParticleTable();

	m_iNumberOfParticlesToBeGenerated = 1;
	m_pParticleDefinition = 0;
	G4ThreeVector hZero(0., 0., 0.);

	m_hParticleMomentumDirection 	= G4ParticleMomentum(1., 0., 0.);
	m_dParticleEnergy 				= 1.0*MeV;
	m_hParticlePosition 			= hZero;
	m_dParticleTime 				= 0.0;
	m_hParticlePolarization 		= hZero;
	m_dParticleCharge				 = 0.0;

	m_hSourcePosType 	= "Volume";
	m_hShape 			= "NULL";
	m_dHalfz 			= 0.;
	m_dRadius 			= 0.;
	m_hCenterCoords 	= hZero;
	m_bConfine 			= false;
	m_hVolumeNames.clear();

	m_hAngDistType 	= "iso";
	m_dMinTheta 	= 0.;
	m_dMaxTheta 	= pi;
	m_dMinPhi 		= 0.;
	m_dMaxPhi 		= twopi;

	m_hEnergyDisType 	= "Mono";
	m_dMonoEnergy 		= 1*MeV;

	// for external file
	in = 0;
	m_hExtStatus = false;

	m_iVerbosityLevel = 0;

	m_pMessenger = new XuerichParticleSourceMessenger(this);
	m_pNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
}

XuerichParticleSource::~XuerichParticleSource()
{
	delete m_pMessenger;
}


/*
//--------------------------------------------------------------------------------------
// SAMPLE PARTICLES FROM EXTERNAL DATA FILE
void
XuerichParticleSource::SetExtFile(G4String hExtFile)
{
	m_hExtFile = hExtFile;
	m_hExtStatus = true;
	m_hShape = "None";
 	m_hAngDistType = "None";
	m_hEnergyDisType = "None";
	ReadExtFile();
}

G4bool
XuerichParticleSource::ReadExtFile()
{
	const double X0 = 0.0;
	const double Y0 = 829.5; 		// right in front of the collimator 
	const double Z0 = 0.0;	

	// read external file
	// energy(keV) | X | Y | Z | cosx | cosy | cosz
	
	ifstream hIn_ext(m_hExtFile.c_str());
	hIn_ext.seekg(in, std::ios::beg);
  
	if(hIn_ext.fail()) {
		G4cout << "Error: cannot open external file " << m_hExtFile << "!" << G4endl;
		return false;
	}

	if(m_iVerbosityLevel >= 1)
 	G4cout << "Source external file: " << m_hExtFile << G4endl;

  	for(int i=0; i<7; i++)
	hIn_ext >> ext_data[i];
  	
  	in=hIn_ext.tellg();
  
    m_pParticleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

	m_dParticleCharge = m_pParticleDefinition->GetPDGCharge();
  
	m_dParticleEnergy = ext_data[0]*keV; // energy is in keV
	
	m_hParticlePosition.setX(X0+ ext_data[1]);
	m_hParticlePosition.setY(Y0+ ext_data[3]);
	m_hParticlePosition.setZ(Z0+ ext_data[2]);
	m_hParticleMomentumDirection.setX(ext_data[4]);
	m_hParticleMomentumDirection.setY(-ext_data[6]);
	m_hParticleMomentumDirection.setZ(ext_data[5]);

	// cout values
	G4cout <<"E  = "<< ext_data[0] << G4endl;
	G4cout <<"X  = "<< ext_data[1] << G4endl;
	G4cout <<"Y  = "<< ext_data[2] << G4endl;
	G4cout <<"Z  = "<< ext_data[3] << G4endl;
	G4cout <<"pX = "<< ext_data[4] << G4endl;
	G4cout <<"pY = "<< ext_data[5] << G4endl;
	G4cout <<"pZ = "<< ext_data[6] << G4endl;

	return true;
}
//--------------------------------------------------------------------------------------
*/


//--------------------------------------------------------------------------------------
// SAMPLE PARTICLES FROM EXTERNAL DATA FILE, UPDATE
void
XuerichParticleSource::SetExtFile(G4String hExtFile)
{
	m_hExtFile = hExtFile;
	m_hExtStatus = true;
	m_hShape = "None";
 	m_hAngDistType = "None";
	m_hEnergyDisType = "None";
	ReadExtFile();
	//first = true;
	in =0;
}

G4bool
XuerichParticleSource::ReadExtFile()
{

	const double X0 = 0.0;
	// Y offset due to collimator thickness is 50.5 mm
	//const double Y0 = 829.5; 	// right in front of the collimator, 88cm 
	//const double Y0 = 649.5; 	// right in front of the collimator, 70cm 
	const double Y0 = 279.5; 	// right in front of the collimator, 28cm 
	const double Z0 = 0.0;	

	// read external file
	// energy(keV) | X | Y | Z | cosx | cosy | cosz
	
	ifstream hIn_ext;
	
		hIn_ext.open(m_hExtFile.c_str());
		hIn_ext.seekg(in, std::ios::beg);

	if(hIn_ext.fail()) {
		G4cout << "Failing to read the file. End of file " << m_hExtFile << ". REWIND" << G4endl;
		in = 0;
		hIn_ext.close();
		hIn_ext.open(m_hExtFile.c_str());
		hIn_ext.clear();
		hIn_ext.seekg(in, std::ios::beg);
	}


/*	if(hIn_ext.fail()) {
		G4cout << "Error: cannot open external file " << m_hExtFile << "!" << G4endl;
		return false;
	}
*/


	if(m_iVerbosityLevel >= 1)
 	G4cout << "Source external file: " << m_hExtFile << G4endl;

	for(int i=0; i<7; i++)
	hIn_ext >> ext_data[i];

  	in=hIn_ext.tellg();
  
  	if(in!=-1){
    m_pParticleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

	m_dParticleCharge = m_pParticleDefinition->GetPDGCharge();
  
	m_dParticleEnergy = ext_data[0]*keV; // energy is in keV
	
	m_hParticlePosition.setX(X0+ ext_data[1]);
	m_hParticlePosition.setY(Y0+ ext_data[3]);
	m_hParticlePosition.setZ(Z0+ ext_data[2]);
	m_hParticleMomentumDirection.setX(ext_data[4]);
	m_hParticleMomentumDirection.setY(-ext_data[6]);
	m_hParticleMomentumDirection.setZ(ext_data[5]);

	// cout values
/*	G4cout <<"E  = "<< ext_data[0] << G4endl;
	G4cout <<"X  = "<< ext_data[1] << G4endl;
	G4cout <<"Y  = "<< ext_data[2] << G4endl;
	G4cout <<"Z  = "<< ext_data[3] << G4endl;
	G4cout <<"pX = "<< ext_data[4] << G4endl;
	G4cout <<"pY = "<< ext_data[5] << G4endl;
	G4cout <<"pZ = "<< ext_data[6] << G4endl;
*/	}
	return true;
}

//--------------------------------------------------------------------------------------




void
XuerichParticleSource::ConfineSourceToVolume(G4String hVolumeList)
{
	stringstream hStream;
	hStream.str(hVolumeList);
	G4String hVolumeName;

	// store all the volume names
	while(!hStream.eof())
	{
		hStream >> hVolumeName;
		m_hVolumeNames.insert(hVolumeName);
	}

	// checks if the selected volumes exist and store all volumes that match
	G4PhysicalVolumeStore *PVStore = G4PhysicalVolumeStore::GetInstance();
	G4bool bFoundAll = true;

	set<G4String> hActualVolumeNames;
	for(set<G4String>::iterator pIt = m_hVolumeNames.begin(); pIt != m_hVolumeNames.end(); pIt++)
	{
		G4String hRequiredVolumeName = *pIt;
		G4bool bMatch = false;

		if(bMatch = (hRequiredVolumeName.last('*') != std::string::npos))
			hRequiredVolumeName = hRequiredVolumeName.strip(G4String::trailing, '*');

		G4bool bFoundOne = false;
		for(G4int iIndex = 0; iIndex < (G4int) PVStore->size(); iIndex++)
		{
			G4String hName = (*PVStore)[iIndex]->GetName();

			if((bMatch && (hName.substr(0, hRequiredVolumeName.size())) == hRequiredVolumeName) || hName == hRequiredVolumeName)
			{
				hActualVolumeNames.insert(hName);
				bFoundOne = true;
			}
		}

		bFoundAll = bFoundAll && bFoundOne;
	}

	if(bFoundAll)
	{
		m_hVolumeNames = hActualVolumeNames;
		m_bConfine = true;

		if(m_iVerbosityLevel >= 1)
			G4cout << "Source confined to volumes: " << hVolumeList << G4endl;

		if(m_iVerbosityLevel >= 2)
		{
			G4cout << "Volume list: " << G4endl;

			for(set<G4String>::iterator pIt = m_hVolumeNames.begin(); pIt != m_hVolumeNames.end(); pIt++)
				G4cout << *pIt << G4endl;
		}
	}
	else if(m_hVolumeNames.empty())
		m_bConfine = false;
	else
	{
		G4cout << " **** Error: One or more volumes do not exist **** " << G4endl;
		G4cout << " Ignoring confine condition" << G4endl;
		m_hVolumeNames.clear();
		m_bConfine = false;
	}
}

void
XuerichParticleSource::GeneratePointSource()
{
	// Generates Points given the point source.
	if(m_hSourcePosType == "Point")
		m_hParticlePosition = m_hCenterCoords;
	else if(m_iVerbosityLevel >= 1)
		G4cout << "Error SourcePosType is not set to Point" << G4endl;
}

void
XuerichParticleSource::GeneratePointsInVolume()
{
	G4ThreeVector RandPos;
	G4double x = 0., y = 0., z = 0.;

	if(m_hSourcePosType != "Volume" && m_iVerbosityLevel >= 1)
		G4cout << "Error SourcePosType not Volume" << G4endl;

	if(m_hShape == "Sphere")
	{
		x = m_dRadius * 2.;
		y = m_dRadius * 2.;
		z = m_dRadius * 2.;
		while(((x * x) + (y * y) + (z * z)) > (m_dRadius * m_dRadius))
		{
			x = G4UniformRand();
			y = G4UniformRand();
			z = G4UniformRand();

			x = (x * 2. * m_dRadius) - m_dRadius;
			y = (y * 2. * m_dRadius) - m_dRadius;
			z = (z * 2. * m_dRadius) - m_dRadius;
		}
	}

	else if(m_hShape == "Cylinder")
	{
		x = m_dRadius * 2.;
		y = m_dRadius * 2.;
		while(((x * x) + (y * y)) > (m_dRadius * m_dRadius))
		{
			x = G4UniformRand();
			y = G4UniformRand();
			z = G4UniformRand();
			x = (x * 2. * m_dRadius) - m_dRadius;
			y = (y * 2. * m_dRadius) - m_dRadius;
			z = (z * 2. * m_dHalfz) - m_dHalfz;
		}
	}

	else if(m_hShape == "None")
		return;

	else
		G4cout << "Error: Volume Shape Does Not Exist" << G4endl;

	RandPos.setX(x);
	RandPos.setY(y);
	RandPos.setZ(z);
	m_hParticlePosition = m_hCenterCoords + RandPos;

}

G4bool
XuerichParticleSource::IsSourceConfined()
{
	// Method to check point is within the volume specified
	if(m_bConfine == false)
		G4cout << "Error: Confine is false" << G4endl;
	G4ThreeVector null(0., 0., 0.);
	G4ThreeVector *ptr;

	ptr = &null;

	// Check m_hParticlePosition is within a volume in our list
	G4VPhysicalVolume *theVolume;

	theVolume = m_pNavigator->LocateGlobalPointAndSetup(m_hParticlePosition, ptr, true);
	G4String theVolName = theVolume->GetName();

	set<G4String>::iterator pIt;
	if((pIt = m_hVolumeNames.find(theVolName)) != m_hVolumeNames.end())
	{
		if(m_iVerbosityLevel >= 1)
			G4cout << "Particle is in volume " << *pIt << G4endl;
		return (true);
	}
	else
		return (false);
}

void
XuerichParticleSource::GenerateIsotropicFlux()
{
	G4double rndm, rndm2;
	G4double px, py, pz;

	G4double sintheta, sinphi, costheta, cosphi;

	rndm = G4UniformRand();
	costheta = std::cos(m_dMinTheta) - rndm * (std::cos(m_dMinTheta) - std::cos(m_dMaxTheta));
	sintheta = std::sqrt(1. - costheta * costheta);

	rndm2 = G4UniformRand();
	m_dPhi = m_dMinPhi + (m_dMaxPhi - m_dMinPhi) * rndm2;
	sinphi = std::sin(m_dPhi);
	cosphi = std::cos(m_dPhi);

	px = -sintheta * cosphi;
	py = -sintheta * sinphi;
	pz = -costheta;

	G4double ResMag = std::sqrt((px * px) + (py * py) + (pz * pz));

	px = px / ResMag;
	py = py / ResMag;
	pz = pz / ResMag;

	m_hParticleMomentumDirection.setX(px);
	m_hParticleMomentumDirection.setY(py);
	m_hParticleMomentumDirection.setZ(pz);

	// m_hParticleMomentumDirection now holds unit momentum vector.
	if(m_iVerbosityLevel >= 2)
		G4cout << "Generating isotropic vector: " <<
			m_hParticleMomentumDirection << G4endl;
}

void
XuerichParticleSource::GenerateMonoEnergetic()
{
	m_dParticleEnergy = m_dMonoEnergy;
}

void
XuerichParticleSource::SetParticleDefinition(G4ParticleDefinition * aParticleDefinition)
{
	m_pParticleDefinition = aParticleDefinition;
	m_dParticleCharge = m_pParticleDefinition->GetPDGCharge();
}

void
XuerichParticleSource::CreatePrimaryParticle(G4PrimaryVertex * vtx)
{
	G4double mass = m_pParticleDefinition->GetPDGMass();
	G4double energy = m_dParticleEnergy + mass;
	G4double pmom = sqrt(energy * energy - mass * mass);
	G4double px = pmom * m_hParticleMomentumDirection.x();
	G4double py = pmom * m_hParticleMomentumDirection.y();
	G4double pz = pmom * m_hParticleMomentumDirection.z();

	primEnergies[NofPrimaries] = m_dParticleEnergy;
	NofPrimaries++;
	emit_energy += m_dParticleEnergy;

	G4PrimaryParticle *particle = new G4PrimaryParticle(m_pParticleDefinition, px, py, pz);
	particle->SetMass(mass);
	particle->SetCharge(m_dParticleCharge);
	particle->SetPolarization(m_hParticlePolarization.x(), m_hParticlePolarization.y(), m_hParticlePolarization.z());
	vtx->SetPrimary(particle);
}


////////////////////////////////////////////////////////////////
// DECAYS
void
XuerichParticleSource::GenerateKr83mDecay(G4PrimaryVertex * vert)
{
	G4double rndm1;
	G4double rndm2;

	rndm1 = G4UniformRand();
	rndm2 = G4UniformRand();
	
	// case 1
	if(rndm1<=0.95)
	{
		if(rndm2<=0.76)
		{	// 1st
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 7.6 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 1.8 * keV;	
			CreatePrimaryParticle(vert);
			// 2nd
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 30.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);
			return;
		}
		else if(rndm2>0.76 && rndm2<=0.85)
		{	// 1st
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 7.6 * keV;	
			CreatePrimaryParticle(vert);
			
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 1.8 * keV;	
			CreatePrimaryParticle(vert);
			// 2nd
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 18.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 10.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);
			return;
		}
		else if(rndm2>0.85)
		{	// 1st
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 7.6 * keV;	
			CreatePrimaryParticle(vert);
			
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 1.8 * keV;	
			CreatePrimaryParticle(vert);
			// 2nd
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 18.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 12.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);
			return;
		}
	}
	// case 2
	else if(rndm1>0.95)
	{
		if(rndm2<=0.76)
		{	// 1st
			particle = particleTable->FindParticle("gamma");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 9.4 * keV;	
			CreatePrimaryParticle(vert);
			// 2nd
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 30.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);
			return;
		}
		else if(rndm2>0.76 && rndm2<=0.85)
		{	// 1st
			particle = particleTable->FindParticle("gamma");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 9.4 * keV;	
			CreatePrimaryParticle(vert);
			// 2nd
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 18.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 10.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);
			return;
		}
		else if(rndm2>0.85)
		{	// 1st
			particle = particleTable->FindParticle("gamma");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 9.4 * keV;	
			CreatePrimaryParticle(vert);
			// 2nd
			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 18.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 12.0 * keV;	
			CreatePrimaryParticle(vert);

			particle = particleTable->FindParticle("e-");
			SetParticleDefinition(particle);
			m_dParticleEnergy = 2.0 * keV;	
			CreatePrimaryParticle(vert);
			return;
		}
	}
}
// DECAYS
////////////////////////////////////////////////////////////




void
XuerichParticleSource::GeneratePrimaryVertex(G4Event * evt)
{

	if(m_pParticleDefinition == 0)
	{
		G4cout << "No particle has been defined!" << G4endl;
		return;
	}

	// Position
	G4bool srcconf = false;
	G4int LoopCount = 0;

	while(srcconf == false)
	{
		if(m_hSourcePosType == "Point")
			GeneratePointSource();
		else if(m_hSourcePosType == "Volume")
			GeneratePointsInVolume();
		else
		{
			G4cout << "Error: SourcePosType undefined" << G4endl;
			G4cout << "Generating point source" << G4endl;
			GeneratePointSource();
		}

		if(m_bConfine == true)
		{
			srcconf = IsSourceConfined();
			// if source in confined srcconf = true terminating the loop
			// if source isnt confined srcconf = false and loop continues
		}
		else if(m_bConfine == false)
			srcconf = true;		// terminate loop

		LoopCount++;
		if(LoopCount == 1000000)
		{
			G4cout << "*************************************" << G4endl;
			G4cout << "LoopCount = 1000000" << G4endl;
			G4cout << "Either the source distribution >> confinement" << G4endl;
			G4cout << "or any confining volume may not overlap with" << G4endl;
			G4cout << "the source distribution or any confining volumes" << G4endl;
			G4cout << "may not exist" << G4endl;
			G4cout << "If you have set confine then this will be ignored" << G4endl;
			G4cout << "for this event." << G4endl;
			G4cout << "*************************************" << G4endl;
			srcconf = true;		//Avoids an infinite loop
		}
	}

	// Angular stuff
	if(m_hAngDistType == "iso")
		GenerateIsotropicFlux();
	else if(m_hAngDistType == "direction")
		SetParticleMomentumDirection(m_hParticleMomentumDirection);
	else if(m_hAngDistType == "None")
		SetParticleMomentumDirection(m_hParticleMomentumDirection);
	else
		G4cout << "Error: AngDistType has unusual value" << G4endl;

	// create a new vertex
	G4PrimaryVertex *vertex = new G4PrimaryVertex(m_hParticlePosition, m_dParticleTime);

	// Energy stuff
	if(m_hEnergyDisType == "Mono")
		GenerateMonoEnergetic();
	else if(m_hEnergyDisType == "Kr83m")
		GenerateKr83mDecay(vertex);
	else if(m_hEnergyDisType != "None")
		G4cout << "Error: EnergyDisType has unusual value" << G4endl;

	if(m_iVerbosityLevel >= 2)
		G4cout << "Creating primaries and assigning to vertex" << G4endl;
	// create new primaries and set them to the vertex
	G4double mass = m_pParticleDefinition->GetPDGMass();
	G4double energy = m_dParticleEnergy + mass;
	G4double pmom = std::sqrt(energy * energy - mass * mass);
	G4double px = pmom * m_hParticleMomentumDirection.x();
	G4double py = pmom * m_hParticleMomentumDirection.y();
	G4double pz = pmom * m_hParticleMomentumDirection.z();

	if(m_iVerbosityLevel >= 1)
	{
		G4cout << "Particle name: " << m_pParticleDefinition->GetParticleName() << G4endl;
		G4cout << "       Energy: " << m_dParticleEnergy << G4endl;
		G4cout << "     Position: " << m_hParticlePosition << G4endl;
		G4cout << "    Direction: " << m_hParticleMomentumDirection << G4endl;
		G4cout << " NumberOfParticlesToBeGenerated: " << m_iNumberOfParticlesToBeGenerated << G4endl;
	}

	for(G4int i = 0; i < m_iNumberOfParticlesToBeGenerated; i++)
	{
		G4PrimaryParticle *particle = new G4PrimaryParticle(m_pParticleDefinition, px, py, pz);
		particle->SetMass(mass);
		particle->SetCharge(m_dParticleCharge);
		particle->SetPolarization(m_hParticlePolarization.x(), m_hParticlePolarization.y(), m_hParticlePolarization.z());
		vertex->SetPrimary(particle);
	}
	evt->AddPrimaryVertex(vertex);
	if(m_iVerbosityLevel > 1)
		G4cout << " Primary Vetex generated " << G4endl;
}

void
XuerichParticleSource::GeneratePrimaryVertexFromTrack(G4Track *pTrack, G4Event *pEvent)
{
	G4double dPX = pTrack->GetMomentum().x();
	G4double dPY = pTrack->GetMomentum().y();
	G4double dPZ = pTrack->GetMomentum().z();

	G4PrimaryVertex *pVertex = new G4PrimaryVertex(pTrack->GetPosition(), m_dParticleTime);

	G4PrimaryParticle *pPrimary = new G4PrimaryParticle(pTrack->GetDefinition(), dPX, dPY, dPZ);
	pPrimary->SetMass(pTrack->GetDefinition()->GetPDGMass());
	pPrimary->SetCharge(pTrack->GetDefinition()->GetPDGCharge());

	pVertex->SetPrimary(pPrimary);

	pEvent->AddPrimaryVertex(pVertex);
}


	