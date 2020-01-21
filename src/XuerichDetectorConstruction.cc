/////////////////////////////////////////////////////////
//													   //
//   ++++++++++++++++++++++++++++++++++++++++++++++++  //
//   + Alexander Kish for XENON-XUERICH experiment	+  //
//   + UZH, 2008-2011								+  //
// 	 ++++++++++++++++++++++++++++++++++++++++++++++++  //
//													   //
/////////////////////////////////////////////////////////

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Polyhedra.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"

#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4VisAttributes.hh"

#include "G4OpBoundaryProcess.hh"
#include "G4SDManager.hh"

#include "G4Colour.hh"
#include "globals.hh"

#include "G4ios.hh"

#include "vector"
#include "numeric"
#include "sstream"
#include "algorithm"
#include "cmath"
#include "cassert"

using std::vector;
using std::stringstream;
using std::max;

#include "XuerichDetectorConstruction.hh"
#include "XuerichDetectorMessenger.hh"

#include "XuerichLXeSensitiveDetector.hh"
#include "XuerichGXeSensitiveDetector.hh"
#include "XuerichPmtSensitiveDetector.hh"
#include "XuerichNaISensitiveDetector.hh"
#include "XuerichDeadXeSensitiveDetector.hh"
#include "XuerichGXeSensitiveDetector.hh"
#include "XuerichTeflonSensitiveDetector.hh"
#include "XuerichSteelSensitiveDetector.hh"
#include "XuerichAlSensitiveDetector.hh"
#include "XuerichLeadSensitiveDetector.hh"
#include "XuerichRestSensitiveDetector.hh"

#include "TMath.h"

XuerichDetectorConstruction::XuerichDetectorConstruction()
{
	m_pDetectorMessenger = new XuerichDetectorMessenger(this);
}

XuerichDetectorConstruction::~XuerichDetectorConstruction()
{
	delete m_pDetectorMessenger;
}


G4VPhysicalVolume *XuerichDetectorConstruction::Construct()
{
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
G4cout 																			<< G4endl;
G4cout << 	" ============================================================"		<< G4endl;
G4cout <<	"|   XUERICH detector model for compton tagging experiment    |"	<< G4endl;
G4cout <<	"| ------------------------                                   |"	<< G4endl;
G4cout <<	"|                           Alexander Kish, UZH 2010-2011    |"	<< G4endl;
G4cout <<	" ============================================================"		<< G4endl;
G4cout <<	"| "																<< G4endl;
////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////

	G4double	density,	// density
    			a,			// atomic weight (g/mole)
    			z;			// atomic number
	G4String	name,		// name
 				symbol;		// symbol
	G4int		ncomponents,
 				natoms;          
	G4double	temperature,
				pressure;    


	//----------	Define elements		--------------------------------
	G4Element *H  = new G4Element(name= "Hydrogen",   symbol= "H",  z= 1.,  a= 1.008   *g/mole);
	G4Element *B  = new G4Element(name= "Boron",	  symbol= "B",  z= 5.,  a= 10.811  *g/mole);
	G4Element *C  = new G4Element(name= "Carbon",     symbol= "C",  z= 6.,  a= 12.011  *g/mole);
	G4Element *N  = new G4Element(name= "Nitrogen",   symbol= "N",  z= 7.,  a= 14.01   *g/mole);
	G4Element *O  = new G4Element(name= "Oxygen",     symbol= "O",  z= 8.,  a= 16.00   *g/mole);
	G4Element *F  = new G4Element(name= "Fluorine",   symbol= "F",  z= 9.,  a= 18.998  *g/mole);
	G4Element *Na = new G4Element(name= "Sodium",     symbol= "Na", z= 11., a= 22.9897 *g/mole);
	G4Element *Al = new G4Element(name= "Aluminium",  symbol= "Al", z= 13., a= 26.98   *g/mole);
	G4Element *Si = new G4Element(name= "Silicon",    symbol= "Si", z= 14., a= 28.0855 *g/mole);
	G4Element *P  = new G4Element(name= "Phosphorus", symbol= "P",  z= 15., a= 30.9738 *g/mole);
	G4Element *S  = new G4Element(name= "Sulphur",    symbol= "S",  z= 16., a= 32.065  *g/mole);
	G4Element *Ar = new G4Element(name= "Argon",	  symbol= "Ar", z= 18., a= 39.948  *g/mole);
	G4Element *Ti = new G4Element(name= "Titanium",	  symbol= "Ti", z= 22., a= 47.867  *g/mole);
	G4Element *Cr = new G4Element(name= "Chromium",   symbol= "Cr", z= 24., a= 51.9962 *g/mole);
	G4Element *Mn = new G4Element(name= "Manganese",  symbol= "Mn", z= 25., a= 54.9381 *g/mole);
	G4Element *Fe = new G4Element(name= "Iron",       symbol= "Fe", z= 26., a= 55.845  *g/mole);
	G4Element *Co = new G4Element(name= "Cobalt",     symbol= "Co", z= 27., a= 58.9332 *g/mole);
	G4Element *Ni = new G4Element(name= "Nickel",     symbol= "Ni", z= 28., a= 58.6934 *g/mole);
	G4Element *Cu = new G4Element(name= "Copper",     symbol= "Cu", z= 29., a= 63.546  *g/mole);
	G4Element *Mo = new G4Element(name= "Molybdenum", symbol= "Mo", z= 42., a= 95.94   *g/mole);
	G4Element *I  = new G4Element(name= "Iodine",     symbol= "I",  z= 53., a= 126.90447 *g/mole);
	G4Element *Xe = new G4Element(name= "Xenon",      symbol= "Xe", z= 54., a= 131.29  *g/mole);
	G4Element *Pb = new G4Element(name= "Lead",       symbol= "Pb", z= 82., a= 207.2   *g/mole);


	////////////////////////////////////////////////////////////////////////////////////////////
	// Define Materials

   // Air
	G4Material *Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

    // Vacuum
	G4Material *Vacuum 	= new G4Material(name = "Vacuum", z = 1., a = 1. *g/mole, density = 1.e-20 *g/cm3, kStateGas, temperature = 0.1 * kelvin, pressure = 1.e-20 * bar);

	// Concrete
	G4Material *Concrete = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");

    // Lead	
    G4Material *Lead_m 	= new G4Material(name = "Lead_m", density = 11.34 *g/cm3, ncomponents =1);
    Lead_m-> AddElement(Pb, 1); 

    // Polyethylene	
	G4Material *Poly 	= new G4Material(name = "Poly", density = 0.95 *g/cm3, ncomponents =2);
	Poly->AddElement(C, 1);
	Poly->AddElement(H, 2);

	// paraffin
	G4Material *Paraffin = G4NistManager::Instance()->FindOrBuildMaterial("G4_PARAFFIN");
	//G4Material *paraffin = new G4Material("paraffin", density=0.93*g/cm3, ncomponents = 2);
	//paraffin->AddElement(H, 14.8605* perCent);
	//paraffin->AddElement(C, 85.1395* perCent);

    // Borated polyethylene (PE 1000 X-Protect, 15% Bortrioxid)
	G4Material *BPE = new G4Material(name = "BPE", density = 1.010 *g/cm3, ncomponents =4);
	BPE->AddElement(C, 0.30);
	BPE->AddElement(H, 0.55);
	BPE->AddElement(B, 0.06);
	BPE->AddElement(O, 0.09);
	
	// wood
	G4Material *wood = new G4Material("wood", density=0.7*g/cm3, ncomponents=3 );
	wood->AddElement(C, 28.57* perCent);
	wood->AddElement(H, 47.61* perCent);
	wood->AddElement(O, 23.82* perCent);
 
	// copper
	G4Material *Copper	= new G4Material(name= "Copper", z= 29., a= 63.5463 *g/mole, 8.92 *g/cm3);
   
  	//----------	Xenon	----------------------------------------
	//G4Material* LXe    = new G4Material(name = "LXe",   2.85*g/cm3, 1, kStateLiquid, 177.15 *kelvin, 1.87 *atmosphere);
	//G4Material* LiqXe  = new G4Material(name = "LiqXe", 2.85*g/cm3, 1, kStateLiquid, 177.15 *kelvin, 1.87 *atmosphere);
	//G4Material* GasXe  = new G4Material("GasXe", 0.005887 * g / cm3, 1, kStateGas, 177.15 * kelvin, 1.5 * atmosphere);
	//LiqXe->AddElement(Xe, 1);	//Liquid Xenon with no scintillation
	//LXe->AddElement(Xe, 1);		//Liquid Xenon with scintillation
	//GasXe->AddElement(Xe, 1);	//Gas Xenon with no Scintillation
	
	// liquid xenon
	G4Material *LXe = new G4Material("LXe", 2.9172*g/cm3, 1, kStateLiquid, 168.15*kelvin, 1.8*atmosphere);
	LXe->AddElement(Xe, 1);

	const G4int iNbEntries = 3;

	G4double pdLXePhotonMomentum[]   	= {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdLXeScintillation[]    	= {0.1,     1.0,     0.1};
	G4double pdLXeRefractiveIndex[]		= {1.63,    1.61,    1.58};
	G4double pdLXeAbsorbtionLength[] 	= {100.*cm, 100.*cm, 100.*cm};
	G4double pdLXeScatteringLength[] 	= {30.*cm,  30.*cm,  30.*cm};
	
	G4MaterialPropertiesTable *pLXePropertiesTable = new G4MaterialPropertiesTable();
	
	pLXePropertiesTable->AddProperty("FASTCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
	pLXePropertiesTable->AddProperty("SLOWCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
	pLXePropertiesTable->AddProperty("RINDEX", pdLXePhotonMomentum, pdLXeRefractiveIndex, iNbEntries);
	pLXePropertiesTable->AddProperty("ABSLENGTH", pdLXePhotonMomentum, pdLXeAbsorbtionLength, iNbEntries);
	pLXePropertiesTable->AddProperty("RAYLEIGH", pdLXePhotonMomentum, pdLXeScatteringLength, iNbEntries);
	
	pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 0./(21.6*eV));
	pLXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 0);
	pLXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
	pLXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
	pLXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);
	
	LXe->SetMaterialPropertiesTable(pLXePropertiesTable);

	// gaseous xenon
	G4Material *GXe = new G4Material("GXe", 0.005887*g/cm3, 1, kStateGas, 173.15*kelvin, 1.5*atmosphere);
	GXe->AddElement(Xe, 1);

	G4double pdGXePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdGXeScintillation[]    = {0.1,     1.0,     0.1};
	G4double pdGXeRefractiveIndex[]  = {1.00,    1.00,    1.00};
	G4double pdGXeAbsorbtionLength[] = {100*m,   100*m,   100*m};
	G4double pdGXeScatteringLength[] = {100*m,   100*m,   100*m};

	G4MaterialPropertiesTable *pGXePropertiesTable = new G4MaterialPropertiesTable();

	pGXePropertiesTable->AddProperty("FASTCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
	pGXePropertiesTable->AddProperty("SLOWCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
	pGXePropertiesTable->AddProperty("RINDEX", pdGXePhotonMomentum, pdGXeRefractiveIndex, iNbEntries);
	pGXePropertiesTable->AddProperty("ABSLENGTH", pdGXePhotonMomentum, pdGXeAbsorbtionLength, iNbEntries);
	pGXePropertiesTable->AddProperty("RAYLEIGH", pdGXePhotonMomentum, pdGXeScatteringLength, iNbEntries);

	pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 0./(keV));
	pGXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 0);
	pGXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
	pGXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
	pGXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);

    GXe->SetMaterialPropertiesTable(pGXePropertiesTable);	

	// argon
	G4Material *LAr		= new G4Material(name = "LAr",		1.40*g/cm3, 1,	kStateLiquid, 87* kelvin,	1. * atmosphere);
	G4Material *LiqAr	= new G4Material(name = "LiqAr",	1.40*g/cm3, 1,	kStateLiquid, 87* kelvin,	1. * atmosphere);
	G4Material *GasAr	= new G4Material(name = "GasAr",	0.001786 * g / cm3, 1, kStateGas, 150.15 * kelvin,1.5 * atmosphere);
	LiqAr	->AddElement(Ar, 1);	//Liquid Argon with no scintillation
	LAr		->AddElement(Ar, 1);	//Liquid Argon with scintillation
	GasAr	->AddElement(Ar, 1);	//Gas Argon with no Scintillation
 
 	// stainless steel
	G4Material *SSteel  = new G4Material(name = "SSteel", density = 7.7 *g/cm3, ncomponents = 3);
	SSteel->AddElement(C,  0.04);
	SSteel->AddElement(Fe, 0.88);
	SSteel->AddElement(Co, 0.08);

  	// mu-metal for NaI light shield
	// source: http://www.drmsmetals.com/data/electronic/almu.html
 	// alternative source: http://www.aircraftmaterialsuk.com/data/electronic/almu.html
	// official source: http://mumetal.co.uk/
	G4Material *MuMetal = new G4Material(name = "MuMetal", density = 8.7 *g/cm3, ncomponents = 5);
	MuMetal->AddElement(Ni,  0.80);
	MuMetal->AddElement(Mo,  0.05);
	MuMetal->AddElement(Si,  0.005);
	MuMetal->AddElement(Cu,  0.0002);
	MuMetal->AddElement(Fe,  0.1448);

    // PTFE
	G4Material *Teflon	= new G4Material("Teflon", 2.16*g/cm3, 2, kStateSolid); // (density for the plain PTFE, i.e. from Boedeker); previously was 2.20g/cm3
	Teflon->AddElement(C, 0.240183);
	Teflon->AddElement(F, 0.759817);

	G4double pdTeflonPhotonMomentum[]  = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdTeflonRefractiveIndex[] = {1.63,    1.61,    1.58};
	G4double pdTeflonReflectivity[]    = {0.95,    0.95,    0.95};
	G4double pdTeflonSpecularLobe[]    = {0.01,    0.01,    0.01};
	G4double pdTeflonSpecularSpike[]   = {0.01,    0.01,    0.01};
	G4double pdTeflonBackscatter[]     = {0.01,    0.01,    0.01};
	G4double pdTeflonEfficiency[]      = {1.0,     1.0,     1.0};
	G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

	pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, iNbEntries);
	pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, iNbEntries);
	pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularLobe, iNbEntries);
	pTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularSpike, iNbEntries);
	pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter, iNbEntries);
	pTeflonPropertiesTable->AddProperty("EFFICIENCY", pdTeflonPhotonMomentum, pdTeflonEfficiency, iNbEntries);
	Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);
		
    // Quartz
	//G4Material* Quartz = new G4Material(name = "Quartz", density = 2.200 *g/cm3, ncomponents = 2);
	//Quartz->AddElement(Si, 1);
	//Quartz->AddElement(O,  2);
	// ref: http://www.sciner.com/Opticsland/FS.htm
	G4Material *Quartz = new G4Material("Quartz", 2.201*g/cm3, 2, kStateSolid, 168.15*kelvin, 1.5*atmosphere);
	Quartz->AddElement(Si, 1);
	Quartz->AddElement(O, 2);

	G4double pdQuartzPhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdQuartzRefractiveIndex[]  = {1.50,    1.56,    1.60};
	G4double pdQuartzAbsorbtionLength[] = {30*m,    30*m,    30*m};

	G4MaterialPropertiesTable *pQuartzPropertiesTable = new G4MaterialPropertiesTable();

	pQuartzPropertiesTable->AddProperty("RINDEX", pdQuartzPhotonMomentum, pdQuartzRefractiveIndex, iNbEntries);
	pQuartzPropertiesTable->AddProperty("ABSLENGTH", pdQuartzPhotonMomentum, pdQuartzAbsorbtionLength, iNbEntries);

	Quartz->SetMaterialPropertiesTable(pQuartzPropertiesTable);

    //---------- Aluminium	------------------------
	G4Material *metalAl	= new G4Material(name = "MetalAluminium", density = 2.700 *g/cm3, ncomponents = 1);
	metalAl->AddElement(Al, 1);
	
	G4Material *Alwire	= new G4Material(name = "Alwire", density = 0.174 *g/cm3, ncomponents = 1);
	Alwire->AddElement(Al, 1);

	G4Material *PhotoCathodeAluminium = new G4Material("PhotoCathodeAluminium", 8.00*g/cm3, 1, kStateSolid);
	PhotoCathodeAluminium->AddElement(Al, 1);


	G4double pdPhotoCathodePhotonMomentum[]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdPhotoCathodeRefractiveIndex[]  = {1.50,    1.56,    1.60};
	G4double pdPhotoCathodeAbsorbtionLength[] = {1.*nm,   1.*nm,   1.*nm};

	G4MaterialPropertiesTable *pPhotoCathodePropertiesTable = new G4MaterialPropertiesTable();

	pPhotoCathodePropertiesTable->AddProperty("RINDEX", pdPhotoCathodePhotonMomentum, pdPhotoCathodeRefractiveIndex, iNbEntries);
	pPhotoCathodePropertiesTable->AddProperty("ABSLENGTH", pdPhotoCathodePhotonMomentum, pdPhotoCathodeAbsorbtionLength, iNbEntries);

	PhotoCathodeAluminium->SetMaterialPropertiesTable(pPhotoCathodePropertiesTable);

    // NaI
    // reference: http://www.optical-components.com/NaI-CsI-crystal.html
	G4Material *NaI 	= new G4Material(name = "NaI", density = 3.67 *g/cm3, ncomponents =2);
	NaI->AddElement(Na, 1);
	NaI->AddElement(I,  1);



////////////////////////////////////////////
//***** Material Optical Properties *****//
	const G4int NUMENTRIES = 3;

//----- Stainless Steel Reflectivity to UV light ---------------------------
	G4double ssteel_PP[NUMENTRIES] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
	G4double ssteel_REFL[NUMENTRIES] = { 0.15, 0.2, 0.15 };
	G4MaterialPropertiesTable *ssteel_mt = new G4MaterialPropertiesTable();

	ssteel_mt->AddProperty("REFLECTIVITY", ssteel_PP, ssteel_REFL, NUMENTRIES);
	SSteel->SetMaterialPropertiesTable(ssteel_mt);

//----- Quartz Window Transparency -----------------------------------------
	G4double quartz_PP[NUMENTRIES]   = { 6.91 * eV, 6.98 * eV, 7.05 * eV };  // lambda range 4 ri
	G4double quartz_RIND[NUMENTRIES] = { 1.50, 1.56, 1.60 };	             // ref index
	G4double quartz_ABSL[NUMENTRIES] = { 30 * m, 30 * m, 30 * m };	         // atten length (remove Quartz)
	G4MaterialPropertiesTable *quartz_mt = new G4MaterialPropertiesTable();

	quartz_mt->AddProperty("RINDEX", quartz_PP, quartz_RIND, NUMENTRIES);
	quartz_mt->AddProperty("ABSLENGTH", quartz_PP, quartz_ABSL, NUMENTRIES);
	Quartz->SetMaterialPropertiesTable(quartz_mt);

//----- PMT Photocathode Absorption ----------------------------------------
	G4double cathmetal_PP[NUMENTRIES] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
	G4double cathmetal_RIND[NUMENTRIES] = { 1.50, 1.56, 1.60 };	                   // ref index
	G4double cathmetal_ABSL[NUMENTRIES] = { 1.e-20 * m, 1.e-20 * m, 1.e-20 * m };  // atten length
	G4MaterialPropertiesTable *cathmetal_mt = new G4MaterialPropertiesTable();

	cathmetal_mt->AddProperty("RINDEX", cathmetal_PP, cathmetal_RIND, NUMENTRIES);
	cathmetal_mt->AddProperty("ABSLENGTH", cathmetal_PP, cathmetal_ABSL, NUMENTRIES);
	metalAl->SetMaterialPropertiesTable(cathmetal_mt);
	
//----- Mesh Transparency --------------------------------------------------
	G4double mesh_PP[NUMENTRIES]   = { 6.91 * eV, 6.98 * eV, 7.05 * eV };	// lambda range 4 ri
	G4double mesh_RIND[NUMENTRIES] = { 1.63, 1.61, 1.58 };	                // ref index, set same as LXe
	G4double absl_len = 1.61615;	                       //atten length, correspond to 94% transparency
	G4double mesh_ABSL[NUMENTRIES] = { absl_len, absl_len, absl_len };
	G4MaterialPropertiesTable *mesh_mt = new G4MaterialPropertiesTable();

	mesh_mt->AddProperty("RINDEX", mesh_PP, mesh_RIND, NUMENTRIES);
	mesh_mt->AddProperty("ABSLENGTH", mesh_PP, mesh_ABSL, NUMENTRIES);
	Alwire->SetMaterialPropertiesTable(mesh_mt);

//////////////////////////////////////////////////////////////////////////////


//--------------	colours		--------------------------------
	G4Colour white  (1.0,	1.0,	1.0);
	G4Colour grey   (0.5,	0.5,	0.5);
	G4Colour lgrey  (.85,	.85,	.85);
	G4Colour red    (1.0,	0.0,	0.0);
	G4Colour lred   (0.75,	0.0,	0.0);
	G4Colour xlred  (0.5,	0.0,	0.0);
	G4Colour cyan   (0.0,	1.0,	1.0);
	G4Colour blue   (0.0,	0.0,	1.0);
	G4Colour lblue  (.5,	0.5,	1.,		1.);
	G4Colour xlblue (.5,	0.5,	1.,		0.2);
	G4Colour magenta(1.0,	0.0,	1.0);
	G4Colour yellow (1.0,	1.0,	0.0);
	G4Colour green  (0.,	.1,		0.);
	G4Colour lgreen (0.0,	.75,	0.0);
	G4Colour xlgreen(0.0,	0.5,	0.0);
	G4Colour brown  (0.7,	0.4,	0.1);
	G4Colour orange (1.0,	0.5,	0.0);

	// rotations
	G4RotationMatrix ZeroRot;
	ZeroRot.rotateX(0. *deg);
	ZeroRot.rotateY(0. *deg);
	ZeroRot.rotateZ(0. *deg);

	G4RotationMatrix RotationXPlus90;
	RotationXPlus90.rotateX(90.*deg);
	RotationXPlus90.rotateY(0.*deg);
	RotationXPlus90.rotateZ(0.*deg);

	G4RotationMatrix RotationXMinus90;
	RotationXMinus90.rotateX(-90.*deg);
	RotationXMinus90.rotateY(0.*deg);
	RotationXMinus90.rotateZ(0.*deg);

	G4double opendeg	= 0.0 *deg;
	G4double closedeg	= 360.0 *deg;

	// Laboratory
	G4double Lab_HalfX = 20.0 *m;
	G4double Lab_HalfY = 20.0 *m;
	G4double Lab_HalfZ = 10.0 *m;

	G4Box *Laboratory = new G4Box("Laboratory", Lab_HalfX, Lab_HalfY, Lab_HalfZ);


	////////////////////////////////////////////////////////////////////////////
	// STAND FOR XUERICH	
	G4double XuStand_thick 				= 6 *mm;	
	G4double XuStand_profile_out 		= 5 *cm;
	G4double XuStand_profile_in 		= XuStand_profile_out - 2*XuStand_thick;

	G4double XuStand_VertBar_height 	= 93 *cm;
	//G4double XuStand_GapTopMid 		= 10 *cm;
	G4double XuStand_GapTopMid 			= 9.15 *cm;

	G4double XuStand_width 	= 50 *cm;
	G4double XuStand_depth 	= 50 *cm;

	G4double XuStand_BottomPlate_thick 	= 1 *cm;
	G4double XuStand_BottomPlate_width 	= 61 *cm;

	G4double XuStand_HorzBarShort_length 	= XuStand_width - 2*XuStand_profile_out;
	
	G4double XuStand_Wing_thick 	= 1.3 *cm;
	G4double XuStand_Wing_length 	= 17 *cm;
	G4double XuStand_Wing_width 	= 5 *cm;

	G4double XuStand_CradlePanel_thick 			= 10 *mm;
	G4double XuStand_CradlePanelLong_length 	= 409 *mm;
	G4double XuStand_CradlePanelLong_width 		= 25 *mm;
	G4double XuStand_CradlePanelShort_length 	= 309 *mm;
	G4double XuStand_CradlePanelShort_width 	= 25 *mm;

	G4double XuStand_CradleBar_height 			= 72.5 *mm;
	G4double XuStand_CradleBarTop_height 		= 62.5 *mm;
	G4double XuStand_CradleBarTop_iD 			= 0.0 *mm;
	G4double XuStand_CradleBarTop_oD 			= 25 *mm;
	G4double XuStand_CradleBarBot_length 		= 49 *mm;
	G4double XuStand_CradleBarBot_width 		= 25 *mm;
	G4double XuStand_CradleBarBot_thick 		= 10 *mm;

	G4double XuStand_Goniometer_thick 			= 10 *mm;
	G4double XuStand_Goniometer_iD 				= 0.0 *mm;
	G4double XuStand_Goniometer_oD 				= 452 *mm;

	G4double XuStand_Cylinder_height 			= 25.0 *cm;
	G4double XuStand_Cylinder_iD 				= 0.00 *mm;
	G4double XuStand_Cylinder_oD 				= 25.0 *mm;
	
	G4Box *XuStand_VertBar_out 	= new G4Box("XuStand_VertBar_out", 	XuStand_profile_out/2, XuStand_profile_out/2, XuStand_VertBar_height/2);
	G4Box *XuStand_VertBar_in 	= new G4Box("XuStand_VertBar_in", XuStand_profile_in/2, XuStand_profile_in/2, XuStand_VertBar_height/2);
	G4ThreeVector XuStand_VertBar_V (0.0, 0.0,	0.0);
	G4Transform3D XuStand_VertBar_T (ZeroRot, XuStand_VertBar_V);	
	G4SubtractionSolid *XuStand_VertBar	= new G4SubtractionSolid( "XuStand_VertBar", XuStand_VertBar_out, XuStand_VertBar_in, XuStand_VertBar_T);
	
	G4Box *XuStand_HorzBarLongX_out = new G4Box("XuStand_HorzBarLongX_out", XuStand_BottomPlate_width/2, XuStand_profile_out/2, XuStand_profile_out/2);
	G4Box *XuStand_HorzBarLongX_in 	= new G4Box("XuStand_HorzBarLongX_in", 	XuStand_BottomPlate_width/2, XuStand_profile_in/2, XuStand_profile_in/2);
	G4ThreeVector XuStand_HorzBarLongX_V (0.0, 0.0,	0.0);
	G4Transform3D XuStand_HorzBarLongX_T (ZeroRot, XuStand_HorzBarLongX_V);	
	G4SubtractionSolid *XuStand_HorzBarLongX	= new G4SubtractionSolid( "XuStand_HorzBarLongX", XuStand_HorzBarLongX_out, XuStand_HorzBarLongX_in, XuStand_HorzBarLongX_T);

	G4Box *XuStand_HorzBarLongY_out 		= new G4Box("XuStand_HorzBarLongY_out", XuStand_profile_out/2, XuStand_depth/2, XuStand_profile_out/2);
	G4Box *XuStand_HorzBarLongY_in 			= new G4Box("XuStand_HorzBarLongY_in", 	XuStand_profile_in/2, XuStand_depth/2, XuStand_profile_in/2);
	G4ThreeVector XuStand_HorzBarLongY_V (0.0, 0.0,	0.0);
	G4Transform3D XuStand_HorzBarLongY_T (ZeroRot, XuStand_HorzBarLongY_V);	
	G4SubtractionSolid *XuStand_HorzBarLongY	= new G4SubtractionSolid( "XuStand_HorzBarLongX", XuStand_HorzBarLongY_out, XuStand_HorzBarLongY_in, XuStand_HorzBarLongY_T);

	G4Box *XuStand_HorzBarShortX_out 		= new G4Box("XuStand_HorzBarShortX_out", XuStand_HorzBarShort_length/2, XuStand_profile_out/2, XuStand_profile_out/2);
	G4Box *XuStand_HorzBarShortX_in 		= new G4Box("XuStand_HorzBarShortX_in", XuStand_HorzBarShort_length/2, XuStand_profile_in/2, XuStand_profile_in/2);
	G4ThreeVector XuStand_HorzBarShortX_V (0.0, 0.0,	0.0);
	G4Transform3D XuStand_HorzBarShortX_T (ZeroRot, XuStand_HorzBarShortX_V);	
	G4SubtractionSolid *XuStand_HorzBarShortX	= new G4SubtractionSolid( "XuStand_HorzBarShortX", XuStand_HorzBarShortX_out, XuStand_HorzBarShortX_in, XuStand_HorzBarShortX_T);

	G4Box *XuStand_HorzBarShortY_out 		= new G4Box("XuStand_HorzBarShortY_out", XuStand_profile_out/2, XuStand_HorzBarShort_length/2, XuStand_profile_out/2);
	G4Box *XuStand_HorzBarShortY_in 		= new G4Box("XuStand_HorzBarShortY_in", XuStand_profile_in/2, XuStand_HorzBarShort_length/2, XuStand_profile_in/2);
	G4ThreeVector XuStand_HorzBarShortY_V (0.0, 0.0,	0.0);
	G4Transform3D XuStand_HorzBarShortY_T (ZeroRot, XuStand_HorzBarShortY_V);
	G4SubtractionSolid *XuStand_HorzBarShortY	= new G4SubtractionSolid( "XuStand_HorzBarShortY", XuStand_HorzBarShortY_out, XuStand_HorzBarShortY_in, XuStand_HorzBarShortY_T);

	G4Box *XuStand_BottomPlate 				= new G4Box("XuStand_HorzBarLong_out", 	XuStand_BottomPlate_width/2, XuStand_BottomPlate_width/2, XuStand_BottomPlate_thick/2);

	G4Box *XuStand_Wing 					= new G4Box("XuStand_Wing", 	XuStand_Wing_width/2, XuStand_Wing_length/2, XuStand_Wing_thick/2);

	G4Tubs *XuStand_Cylinder				= new G4Tubs("XuStand_Cylinder", XuStand_Cylinder_iD/2, XuStand_Cylinder_oD/2, XuStand_Cylinder_height/2, opendeg, closedeg);
	G4Box  *XuStand_CradlePanelLong 		= new G4Box("XuStand_CradlePanelLong", 	XuStand_CradlePanelLong_width/2, XuStand_CradlePanelLong_length/2, XuStand_CradlePanel_thick/2);
	G4Box  *XuStand_CradlePanelShort		= new G4Box("XuStand_CradlePanelShort", XuStand_CradlePanelShort_width/2, XuStand_CradlePanelShort_length/2, XuStand_CradlePanel_thick/2);

	G4Tubs *XuStand_CradleBarTop			= new G4Tubs("XuStand_CradleBarTop", XuStand_CradleBarTop_iD/2, XuStand_CradleBarTop_oD/2, XuStand_CradleBarTop_height/2, opendeg, closedeg);
	G4Box  *XuStand_CradleBarBot 			= new G4Box("XuStand_CradleBarBot", XuStand_CradleBarBot_width/2, XuStand_CradleBarBot_length/2, XuStand_CradleBarBot_thick/2);
	G4ThreeVector XuStand_CradleBarFront_V	(0.0, XuStand_CradleBarTop_oD/2-XuStand_CradleBarBot_length/2, -XuStand_CradleBarTop_height/2-XuStand_CradleBarBot_thick/2);
	G4Transform3D XuStand_CradleBarFront_T	(ZeroRot, XuStand_CradleBarFront_V);
	G4UnionSolid *XuStand_CradleBarFront	= new G4UnionSolid( "XuStand_CradleBarFront",	XuStand_CradleBarTop,	XuStand_CradleBarBot,	XuStand_CradleBarFront_T);
	G4ThreeVector XuStand_CradleBarBack_V	(0.0, -XuStand_CradleBarTop_oD/2+XuStand_CradleBarBot_length/2, -XuStand_CradleBarTop_height/2-XuStand_CradleBarBot_thick/2);
	G4Transform3D XuStand_CradleBarBack_T	(ZeroRot, XuStand_CradleBarBack_V);
	G4UnionSolid *XuStand_CradleBarBack		= new G4UnionSolid( "XuStand_CradleBarBack",	XuStand_CradleBarTop,	XuStand_CradleBarBot, XuStand_CradleBarBack_T);

	G4Tubs *XuStand_Goniometer				= new G4Tubs("XuStand_Goniometer", XuStand_Goniometer_iD/2, XuStand_Goniometer_oD/2, XuStand_Goniometer_thick/2, opendeg, closedeg);


	/////////////////////////////////
	// XUERICH DETECTOR
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Cryostat: Inner Can
	G4double InnerCanTube_WallThick		= 2.0 *mm;
	G4double InnerCanTube_height		= 180.0 *mm;
	G4double InnerCanTube_iR			= 74.5 *mm;
	G4double InnerCanTube_oR			= InnerCanTube_iR + InnerCanTube_WallThick;
	G4double InnerCanBottom_thick		= InnerCanTube_WallThick;
	G4double InnerCanBottom_iR			= 0.0 *mm;
	G4double InnerCanBottom_oR			= InnerCanTube_oR;	
	G4double InnerCanLowerFlange_iR		= InnerCanTube_iR;
	G4double InnerCanLowerFlange_oR		= 101.5 *mm;
	G4double InnerCanLowerFlange_thick	= 22.0 *mm;
	G4double InnerCanUpperFlange_iR		= 0.0 *mm;
	G4double InnerCanUpperFlange_oR		= InnerCanLowerFlange_oR;
	G4double InnerCanUpperFlange_thick	= InnerCanLowerFlange_thick;
	G4double InnerCan_height			= InnerCanTube_height + InnerCanBottom_thick + InnerCanLowerFlange_thick + InnerCanUpperFlange_thick;
	
	G4Tubs *InnerCanTube 		= new G4Tubs("InnerCanTube", InnerCanTube_iR, InnerCanTube_oR, InnerCanTube_height/2, opendeg, closedeg);
	G4Tubs *InnerCanBottom 		= new G4Tubs("InnerCanBottom", InnerCanBottom_iR, InnerCanBottom_oR, InnerCanBottom_thick/2, opendeg, closedeg);
	G4Tubs *InnerCanLowerFlange	= new G4Tubs("InnerCanLowerFlange", InnerCanLowerFlange_iR, InnerCanLowerFlange_oR, InnerCanLowerFlange_thick/2, opendeg, closedeg);
	G4Tubs *InnerCanUpperFlange	= new G4Tubs("InnerCanUpperFlange", InnerCanUpperFlange_iR, InnerCanUpperFlange_oR, InnerCanUpperFlange_thick/2, opendeg, closedeg);
	
	// Cryostat: Outer Can
	G4double OuterCanTube_WallThick		= 3.0 *mm;
	G4double OuterCanTube_height		= 380.0 *mm;
	G4double OuterCanTube_iR			= 125.0 *mm;
	G4double OuterCanTube_oR			= OuterCanTube_iR + OuterCanTube_WallThick;
	G4double OuterCanBottom_thick		= OuterCanTube_WallThick;
	G4double OuterCanBottom_iR			= 0.0 *mm;
	G4double OuterCanBottom_oR			= OuterCanTube_oR;
	G4double OuterCanLowerFlange_thick	= 12.0 *mm;
	G4double OuterCanLowerFlange_iR		= OuterCanTube_iR;
	G4double OuterCanLowerFlange_oR		= 145.0 *mm;
	G4double OuterCanUpperFlange_thick	= 12.0 *mm;
	G4double OuterCanUpperFlange_iR		= 0.0 *mm;
	G4double OuterCanUpperFlange_oR		= 145.0 *mm;
	G4double OuterCan_height			= OuterCanTube_height + OuterCanBottom_thick + OuterCanLowerFlange_thick + OuterCanUpperFlange_thick;
	
	G4Tubs *OuterCanTube		= new G4Tubs("OuterCanTube", OuterCanTube_iR, OuterCanTube_oR, OuterCanTube_height/2, opendeg, closedeg);
	G4Tubs *OuterCanBottom		= new G4Tubs("OuterCanBottom", OuterCanBottom_iR, OuterCanBottom_oR, OuterCanBottom_thick/2, opendeg, closedeg);
	G4Tubs *OuterCanLowerFlange	= new G4Tubs("OuterCanLowerFlange", OuterCanLowerFlange_iR, OuterCanLowerFlange_oR, OuterCanLowerFlange_thick/2, opendeg, closedeg);
	G4Tubs *OuterCanUpperFlange	= new G4Tubs("OuterCanUpperFlange", OuterCanUpperFlange_iR, OuterCanUpperFlange_oR, OuterCanUpperFlange_thick/2, opendeg, closedeg);

	//--------------	Cryostat: Vacuum	------------------------
	G4double CryostatVacuum_iR		= 0.0 *mm;
	G4double CryostatVacuum_oR		= OuterCanTube_iR;
	G4double CryostatVacuum_height	= OuterCanTube_height;
	G4Tubs *CryostatVacuum	= new G4Tubs("CryostatVacuum", CryostatVacuum_iR, CryostatVacuum_oR, CryostatVacuum_height/2, opendeg, closedeg);

	// Cryostat: Al can
	G4double AlCanTube_WallThick	= 5.5 *mm;
	G4double AlCanTube_height		= 245.0 *mm;
	G4double AlCanTube_oD			= 215.0 *mm;
	G4double AlCanTube_iD			= AlCanTube_oD - 2*AlCanTube_WallThick;
	G4double AlCanBottom_thick		= 8.5 *mm;
	G4double AlCanBottom_iD			= 0.0;
	G4double AlCanBottom_oD			= AlCanTube_oD;
	G4double AlCan_height			= AlCanTube_height + AlCanBottom_thick;
	
	G4Tubs *AlCanTube	= new G4Tubs("AlCanTube", AlCanTube_iD/2, AlCanTube_oD/2, AlCanTube_height/2, opendeg, closedeg);
	G4Tubs *AlCanBottom	= new G4Tubs("AlCanBottom", AlCanBottom_iD/2, AlCanBottom_oD/2, AlCanBottom_thick/2, opendeg, closedeg);
	
	// Copper Finger
	G4double CopperFinger_iD 		= 0.0 *mm;
	G4double CopperFinger_oD 		= 25 *mm;
	G4double CopperFinger_height	= 525 *mm;	
	G4Tubs *CopperFinger = new G4Tubs("CopperFinger", CopperFinger_iD/2, CopperFinger_oD/2, CopperFinger_height/2, opendeg, closedeg);
	
	// dewar
	//G4double Dewar_height 		= 73.0 *cm;
	G4double Dewar_height 			= 69.0 *cm;
	G4double Dewar_thick 			= 1.5 *mm;
	G4double DewarVacuum_thick 		= 1.0 *cm;
	G4double DewarOutTube_height 	= Dewar_height - 2*Dewar_thick;
	G4double DewarInTube_height 	= Dewar_height - 3*Dewar_thick - DewarVacuum_thick;
	G4double DewarOutTube_oD 		= 25.0 *cm;
	G4double DewarOutTube_iD 		= DewarOutTube_oD - 2*Dewar_thick;
	G4double DewarOutBot_oD 		= DewarOutTube_oD;
	G4double DewarOutBot_iD 		= 0.0;
	G4double DewarInTube_oD 		= DewarOutTube_iD - 2*DewarVacuum_thick;
	G4double DewarInTube_iD 		= DewarInTube_oD - 2*Dewar_thick;
	G4double DewarInBot_oD 			= DewarInTube_oD;
	G4double DewarInBot_iD 			= 0.0;
	G4double DewarTop_iD			= DewarInTube_iD;
	G4double DewarTop_oD			= DewarOutTube_oD;
	
	G4Tubs *DewarOutTube 	= new G4Tubs("DewarOutTube", DewarOutTube_iD/2, DewarOutTube_oD/2, DewarOutTube_height/2, opendeg, closedeg);
	G4Tubs *DewarOutBot 	= new G4Tubs("DewarOutBot", DewarOutBot_iD/2, DewarOutBot_oD/2, Dewar_thick/2, opendeg, closedeg);
	G4Tubs *DewarInTube 	= new G4Tubs("DewarInTube", DewarInTube_iD/2, DewarInTube_oD/2, DewarInTube_height/2, opendeg, closedeg);
	G4Tubs *DewarInBot 		= new G4Tubs("DewarInBot", DewarInBot_iD/2, DewarInBot_oD/2, Dewar_thick/2, opendeg, closedeg);
	G4Tubs *DewarTop 		= new G4Tubs("DewarTop", DewarTop_iD/2, DewarTop_oD/2, Dewar_thick/2, opendeg, closedeg);
		
	// PIPES, FLANGES
	G4double InnerPipes_height	= 101 *mm;
	G4double Pipes_thick		= 2.0 *mm;

	// pipe 1 (8), pressure gauge
	G4double Pipe1_oD 			= 38 *mm;
	G4double Pipe1_iD 			= Pipe1_oD - 2*Pipes_thick;
	G4double Pipe1_height 		= 250 *mm;
	G4double Pipe1flange_thick 	= 20 *mm;
	G4double Pipe1flange_oD 	= 70 *mm;
	G4double Pipe1flange_iD 	= Pipe1_oD;
	G4double Pipe1flange_offset = 210 *mm;

	// pipe 2 (5)
	G4double Pipe2_oD 			= 16 *mm;
	G4double Pipe2_iD 			= Pipe2_oD - 2*Pipes_thick;
	G4double Pipe2_height 		= 190 *mm;
	G4double Pipe2flange_thick 	= 14 *mm;
	G4double Pipe2flange_oD 	= 33 *mm;
	G4double Pipe2flange_iD 	= Pipe2_oD;
	G4double Pipe2flange_offset = 100 *mm;

	// pipe 3 (6)
	G4double Pipe3_oD 			= 38 *mm;
	G4double Pipe3_iD 			= Pipe3_oD - 2*Pipes_thick;
	G4double Pipe3_height 		= 180 *mm;
	G4double Pipe3flange_thick 	= 20 *mm;
	G4double Pipe3flange_oD 	= 70 *mm;
	G4double Pipe3flange_iD 	= Pipe3_oD;
	G4double Pipe3flange_offset = 100 *mm;

	// pipe 4 (11)
	G4double Pipe4_oD 			= 38 *mm;
	G4double Pipe4_iD 			= Pipe4_oD - 2*Pipes_thick;
	G4double Pipe4_height 		= 75 *mm;
	G4double Pipe4flange_thick 	= 20 *mm;
	G4double Pipe4flange_oD 	= 70 *mm;
	G4double Pipe4flange_iD 	= Pipe4_oD;
	G4double Pipe4flange_offset = 45 *mm;

	// pipe 5 (2), pumping port
	G4double Pipe5_oD 			= 38 *mm;
	G4double Pipe5_iD 			= Pipe5_oD - 2*Pipes_thick;
	G4double Pipe5_height 		= 270 *mm;
	G4double Pipe5flange_thick 	= 20 *mm;
	G4double Pipe5flange_oD 	= 70 *mm;
	G4double Pipe5flange_iD 	= Pipe5_oD;
	G4double Pipe5flange1_offset= 95 *mm;
	G4double Pipe5flange2_offset= 230 *mm;

	// pipe 6 (1)
	G4double Pipe6_oD 			= 16 *mm;
	G4double Pipe6_iD 			= Pipe6_oD - 2*Pipes_thick;
	G4double Pipe6_height 		= 190 *mm;
	G4double Pipe6flange_thick 	= 14 *mm;
	G4double Pipe6flange_oD 	= 33 *mm;
	G4double Pipe6flange_iD 	= Pipe6_oD;
	G4double Pipe6flange_offset	= 100 *mm;

	// pipe 7 (3), xenon line
	G4double Pipe7_oD 			= 38 *mm;
	G4double Pipe7_iD 			= Pipe7_oD - 2*Pipes_thick;
	G4double Pipe7_height 		= 180 *mm;
	G4double Pipe7flange_thick 	= 20 *mm;
	G4double Pipe7flange_oD 	= 70 *mm;
	G4double Pipe7flange_iD 	= Pipe7_oD;

	// pipe 8 (4), PMT HV
	G4double Pipe8_oD 			= 16 *mm;
	G4double Pipe8_iD 			= Pipe8_oD - 2*Pipes_thick;
	G4double Pipe8_height 		= 190 *mm;
	G4double Pipe8flange_thick 	= 14 *mm;
	G4double Pipe8flange_oD 	= 33 *mm;
	G4double Pipe8flange_iD 	= Pipe8_oD;
	G4double Pipe8flange_offset	= 100 *mm;

	// pipe 9 (7)
	G4double Pipe9_oD 			= 38 *mm;
	G4double Pipe9_iD 			= Pipe9_oD - 2*Pipes_thick;
	G4double Pipe9_height 		= 180 *mm;
	G4double Pipe9flange_thick 	= 20 *mm;
	G4double Pipe9flange_oD 	= 70 *mm;
	G4double Pipe9flange_iD 	= Pipe9_oD;
	G4double Pipe9flange_offset	= 100 *mm;

	// pipe 10 (9)
	G4double Pipe10_oD 			= 38 *mm;
	G4double Pipe10_iD 			= Pipe10_oD - 2*Pipes_thick;
	G4double Pipe10_height 		= 75 *mm;
	G4double Pipe10flange_thick = 20 *mm;
	G4double Pipe10flange_oD 	= 70 *mm;
	G4double Pipe10flange_iD 	= Pipe10_oD;
	G4double Pipe10flange_offset= 45 *mm;

	// pipe 11 (10)
	G4double Pipe11_oD 			= 38 *mm;
	G4double Pipe11_iD 			= Pipe11_oD - 2*Pipes_thick;
	G4double Pipe11_height 		= 75 *mm;
	G4double Pipe11flange_thick = 20 *mm;
	G4double Pipe11flange_oD 	= 70 *mm;
	G4double Pipe11flange_iD 	= Pipe11_oD;
	G4double Pipe11flange_offset= 45 *mm;

	G4Tubs *Pipe1in		= new G4Tubs("Pipe1in", 	Pipe1_iD/2, 	Pipe1_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe2in		= new G4Tubs("Pipe2in", 	Pipe2_iD/2, 	Pipe2_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe3in		= new G4Tubs("Pipe3in", 	Pipe3_iD/2, 	Pipe3_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe4in		= new G4Tubs("Pipe4in", 	Pipe4_iD/2, 	Pipe4_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe5in		= new G4Tubs("Pipe5in", 	Pipe5_iD/2, 	Pipe5_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe6in		= new G4Tubs("Pipe6in", 	Pipe6_iD/2, 	Pipe6_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe7in		= new G4Tubs("Pipe7in", 	Pipe7_iD/2, 	Pipe7_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe8in		= new G4Tubs("Pipe8in", 	Pipe8_iD/2, 	Pipe8_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe9in		= new G4Tubs("Pipe9in", 	Pipe9_iD/2, 	Pipe9_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe10in	= new G4Tubs("Pipe10in", 	Pipe10_iD/2,	Pipe10_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe11in	= new G4Tubs("Pipe11in", 	Pipe11_iD/2, 	Pipe11_oD/2, 	InnerPipes_height/2, 	opendeg, closedeg);

	G4Tubs *Pipe1out	= new G4Tubs("Pipe1out", 	Pipe1_iD/2, 	Pipe1_oD/2, 	Pipe1_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe2out	= new G4Tubs("Pipe2out", 	Pipe2_iD/2, 	Pipe2_oD/2, 	Pipe2_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe3out	= new G4Tubs("Pipe3out", 	Pipe3_iD/2, 	Pipe3_oD/2, 	Pipe3_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe4out	= new G4Tubs("Pipe4out", 	Pipe4_iD/2, 	Pipe4_oD/2, 	Pipe4_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe5out	= new G4Tubs("Pipe5out", 	Pipe5_iD/2, 	Pipe5_oD/2, 	Pipe5_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe6out	= new G4Tubs("Pipe6out", 	Pipe6_iD/2, 	Pipe6_oD/2, 	Pipe6_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe7out	= new G4Tubs("Pipe7out", 	Pipe7_iD/2, 	Pipe7_oD/2, 	Pipe7_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe8out	= new G4Tubs("Pipe8out", 	Pipe8_iD/2, 	Pipe8_oD/2, 	Pipe8_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe9out	= new G4Tubs("Pipe9out", 	Pipe9_iD/2, 	Pipe9_oD/2, 	Pipe9_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe10out	= new G4Tubs("Pipe10out", 	Pipe10_iD/2,	Pipe10_oD/2,	Pipe10_height/2, 	opendeg, closedeg);
	G4Tubs *Pipe11out	= new G4Tubs("Pipe11out", 	Pipe11_iD/2, 	Pipe11_oD/2, 	Pipe11_height/2, 	opendeg, closedeg);

	G4Tubs *Pipe1flange	= new G4Tubs("Pipe1flange", 	Pipe1flange_iD/2, 	Pipe1flange_oD/2, 	Pipe1flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe2flange	= new G4Tubs("Pipe2flange", 	Pipe2flange_iD/2, 	Pipe2flange_oD/2, 	Pipe2flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe3flange	= new G4Tubs("Pipe3flange", 	Pipe3flange_iD/2, 	Pipe3flange_oD/2, 	Pipe3flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe4flange	= new G4Tubs("Pipe4flange", 	Pipe4flange_iD/2, 	Pipe4flange_oD/2, 	Pipe4flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe5flange	= new G4Tubs("Pipe5flange", 	Pipe5flange_iD/2, 	Pipe5flange_oD/2, 	Pipe5flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe6flange	= new G4Tubs("Pipe6flange", 	Pipe6flange_iD/2, 	Pipe6flange_oD/2, 	Pipe6flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe8flange	= new G4Tubs("Pipe8flange", 	Pipe8flange_iD/2, 	Pipe8flange_oD/2, 	Pipe8flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe9flange	= new G4Tubs("Pipe9flange", 	Pipe9flange_iD/2, 	Pipe9flange_oD/2, 	Pipe9flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe10flange= new G4Tubs("Pipe10flange", 	Pipe10flange_iD/2,	Pipe10flange_oD/2,	Pipe10flange_thick/2, 	opendeg, closedeg);
	G4Tubs *Pipe11flange= new G4Tubs("Pipe11flange", 	Pipe11flange_iD/2, 	Pipe11flange_oD/2, 	Pipe11flange_thick/2, 	opendeg, closedeg);
	
	//--------------	SteelHolder			------------------------
	G4double SteelHolder_thick				= 3.0 *mm;
	G4double SteelHolder_height				= 38.0 *mm;
	G4double SteelHolderTube_iR				= 28.5 *mm;	
	G4double SteelHolderTube_oR				= SteelHolderTube_iR + SteelHolder_thick;
	G4double SteelHolderTube_height			= SteelHolder_height - 2*SteelHolder_thick;
	G4double SteelHolderUpperSquare_side	= 65.0 *mm;
	G4double SteelHolderUpperSquare_thick	= SteelHolder_thick;
	G4double SteelHolderLowerSquare_side	= 65.0 *mm;
	G4double SteelHolderLowerSquare_thick	= SteelHolder_thick;
	G4double SteelHolderUpperCutout_iR		= 0.0 *mm;
	G4double SteelHolderUpperCutout_oR		= SteelHolderTube_iR;
	G4double SteelHolderUpperCutout_thick	= SteelHolder_thick;
	G4double SteelHolderLowerCutout_iR		= 0.0 *mm;
	G4double SteelHolderLowerCutout_oR		= SteelHolderTube_iR;
	G4double SteelHolderLowerCutout_thick	= SteelHolder_thick;

	G4Tubs *SteelHolderTube			= new G4Tubs("SteelHolderTube", SteelHolderTube_iR, SteelHolderTube_oR, SteelHolderTube_height/2, opendeg, closedeg);
	G4Box *SteelHolderUpperSquare 	= new G4Box("SteelHolderUpperSquare", SteelHolderUpperSquare_side/2, SteelHolderUpperSquare_side/2, SteelHolderUpperSquare_thick/2);
	G4Box *SteelHolderLowerSquare	= new G4Box("SteelHolderLowerSquare", SteelHolderLowerSquare_side/2, SteelHolderLowerSquare_side/2, SteelHolderLowerSquare_thick/2);
	G4Tubs *SteelHolderUpperCutout	= new G4Tubs("SteelHolderUpperCutout", SteelHolderUpperCutout_iR, SteelHolderUpperCutout_oR, SteelHolderUpperCutout_thick/2, opendeg, closedeg);
	G4Tubs *SteelHolderLowerCutout	= new G4Tubs("SteelHolderLowerCutOut", SteelHolderLowerCutout_iR, SteelHolderLowerCutout_oR, SteelHolderLowerCutout_thick/2, opendeg, closedeg);

	G4ThreeVector SteelHolderUpperCutout_V	(0.0, 0.0,	0.0);	
	G4Transform3D SteelHolderUpperCutout_T	(ZeroRot, 	SteelHolderUpperCutout_V);

	G4ThreeVector SteelHolderLowerCutout_V	(0.0, 0.0,	0.0);	
	G4Transform3D SteelHolderLowerCutout_T	(ZeroRot, 	SteelHolderLowerCutout_V);

	G4SubtractionSolid
	* SteelHolderUpperSquare1
	= new G4SubtractionSolid(	"SteelHolderUpperSquare1",	SteelHolderUpperSquare,	SteelHolderUpperCutout,		SteelHolderUpperCutout_T);

	G4SubtractionSolid
	* SteelHolderLowerSquare1
	= new G4SubtractionSolid(	"SteelHolderLowerSquare1",	SteelHolderLowerSquare,	SteelHolderLowerCutout,		SteelHolderLowerCutout_T);

	//--------------	Supports			------------------------
	G4double Support_height		= 34.0 *mm;
	G4double Support_iR			= 0.0 *mm;
	G4double Support_oR			= 8.0 *mm;
	G4double Washer_iR			= 0.0 *mm;
	G4double Washer_oR			= 6.0 *mm;
	G4double Washer_thick		= 2.0 *mm;
	
	G4Tubs *Support = new G4Tubs("Support", Support_iR, Support_oR, Support_height/2, opendeg, closedeg);
	G4Tubs *Washer	= new G4Tubs("Washer", Washer_iR, Washer_oR, Washer_thick/2, opendeg, closedeg);

	//--------------	Top PMT Clamp		------------------------
	G4double TopClamp_iR	= 22.5 *mm;
	G4double TopClamp_oR	= 30.0 *mm;
	G4double TopClamp_thick	= 2.0 *mm;
	
	G4Tubs *TopClamp = new G4Tubs("TopClamp", TopClamp_iR, TopClamp_oR, TopClamp_thick/2, opendeg, closedeg);
	
	//--------------	Top PMT Holder		------------------------
	G4double TopHolderTube_oR				= 60.0 *mm;
	G4double TopHolderTube_iR				= 0.0 *mm;
	G4double TopHolderTube_height			= 33.0 *mm;
	
	G4double TopHolderUpperTube_iR			= 28.5 *mm;	
	G4double TopHolderUpperTube_oR			= TopHolderTube_oR;
	G4double TopHolderUpperTube_height		= 28.0 *mm;
	G4double TopHolderLowerTube_iR			= 17.5 *mm;
	G4double TopHolderLowerTube_oR			= TopHolderTube_oR;
	G4double TopHolderLowerTube_height		= TopHolderTube_height - TopHolderUpperTube_height;
	
	G4double TopHolderSquare_side			= 58.5 *mm;
	G4double TopHolderSquare_thick			= 2.5 *mm;
	G4double TopHolderSquareCutout_iR		= 0.0 *mm;
	G4double TopHolderSquareCutout_oR		= 17.5 *mm;
	G4double TopHolderSquareCutout_thick	= TopHolderSquare_thick;	
	
	G4double TopHolder_height				= TopHolderTube_height + TopHolderSquare_thick;
	
	G4Tubs *TopHolderUpperTube 		= new G4Tubs("TopHolderUpperTube", TopHolderUpperTube_iR, TopHolderUpperTube_oR, TopHolderUpperTube_height/2, opendeg, closedeg);
	G4Tubs *TopHolderLowerTube		= new G4Tubs("TopHolderLowerTube", TopHolderLowerTube_iR, TopHolderLowerTube_oR, TopHolderLowerTube_height/2, opendeg, closedeg);
	G4Box *TopHolderSquare			= new G4Box("TopHolderSquare", TopHolderSquare_side/2, TopHolderSquare_side/2, TopHolderSquare_thick/2);
	G4Tubs *TopHolderSquareCutout	= new G4Tubs("TopHolderSquareCutout", TopHolderSquareCutout_iR, TopHolderSquareCutout_oR, TopHolderSquareCutout_thick/2, opendeg, closedeg);
	
	G4ThreeVector TopHolderSquareCutout_V	(0.0, 0.0,  0.0);
	G4Transform3D TopHolderSquareCutout_T	(ZeroRot, TopHolderSquareCutout_V);

	G4SubtractionSolid *TopHolderSquare1 = new G4SubtractionSolid("TopHolderSquare1",	TopHolderSquare,	TopHolderSquareCutout,	TopHolderSquareCutout_T);

	//--------------	Extraction Spacer	------------------------
	G4double ExtractionSpacer_iR		= 17.5 *mm;
	G4double ExtractionSpacer_oR		= 60.0 *mm;
	G4double ExtractionSpacer_height	= 5.0 *mm;
	
	G4double ExtractionHalfSpacer_iR		= ExtractionSpacer_iR;
	G4double ExtractionHalfSpacer_oR		= ExtractionSpacer_oR;
	G4double ExtractionHalfSpacer_height	= ExtractionSpacer_height/2;
	
	G4Tubs *ExtractionHalfSpacer = new G4Tubs("ExtractionHalfSpacer", ExtractionHalfSpacer_iR, ExtractionHalfSpacer_oR, ExtractionHalfSpacer_height/2, opendeg, closedeg);

	// Drift Spacer (normal big spacer)
	G4double DriftSpacerTube_iR				= 17.5 *mm;
	G4double DriftSpacerTube_oR				= 60.0 *mm;
	G4double DriftSpacerTube_height			= 24.0 *mm;	

	G4double DriftSpacerUpperTube_iR		= DriftSpacerTube_iR;
	G4double DriftSpacerUpperTube_oR		= DriftSpacerTube_oR;
	G4double DriftSpacerUpperTube_height	= 11.0 *mm;
	
	G4double DriftSpacerMiddleTube_iR		= DriftSpacerTube_iR;
	G4double DriftSpacerMiddleTube_oR		= 30.0 *mm;
	G4double DriftSpacerMiddleTube_height	= 2.0 *mm;
	
	
	G4double DriftSpacerLowerTube_iR		= DriftSpacerTube_iR;
	G4double DriftSpacerLowerTube_oR		= DriftSpacerTube_oR;
	G4double DriftSpacerLowerTube_height	= 11.0 *mm;

	G4double DriftSpacerUpperSquare_side	= 58.5 *mm;
	G4double DriftSpacerUpperSquare_thick	= 2.5 *mm;
	G4double DriftSpacerLowerSquare_side	= 58.5 *mm;
	G4double DriftSpacerLowerSquare_thick	= 2.5 *mm;

	G4double DriftSpacerUpperSquareCutout_iR	= 0.0 *mm;
	G4double DriftSpacerUpperSquareCutout_oR	= 17.5 *mm;
	G4double DriftSpacerUpperSquareCutout_thick	= DriftSpacerUpperSquare_thick;
	G4double DriftSpacerLowerSquareCutout_iR	= 0.0 *mm;
	G4double DriftSpacerLowerSquareCutout_oR	= 17.5 *mm;
	G4double DriftSpacerLowerSquareCutout_thick	= DriftSpacerLowerSquare_thick;
	
	G4Tubs *DriftSpacerUpperTube 		= new G4Tubs("DriftSpacerUpperTube",	DriftSpacerUpperTube_iR, DriftSpacerUpperTube_oR, DriftSpacerUpperTube_height/2, opendeg, closedeg);
	G4Tubs *DriftSpacerMiddleTube		= new G4Tubs("DriftSpacerMiddleTube",	DriftSpacerMiddleTube_iR, DriftSpacerMiddleTube_oR, DriftSpacerMiddleTube_height/2, opendeg, closedeg);
	G4Tubs *DriftSpacerLowerTube		= new G4Tubs("DriftSpacerLowerTube",	DriftSpacerLowerTube_iR, DriftSpacerLowerTube_oR, DriftSpacerLowerTube_height/2, opendeg, closedeg);
	G4Box *DriftSpacerUpperSquare		= new G4Box("DriftSpacerUpperSquare",	DriftSpacerUpperSquare_side/2, DriftSpacerUpperSquare_side/2, DriftSpacerUpperSquare_thick/2);
	G4Box *DriftSpacerLowerSquare		= new G4Box("DriftSpacerLowerSquare",	DriftSpacerLowerSquare_side/2, DriftSpacerLowerSquare_side/2, DriftSpacerLowerSquare_thick/2);
	G4Tubs *DriftSpacerUpperSquareCutout= new G4Tubs("DriftSpacerUpperSquareCutout",	DriftSpacerUpperSquareCutout_iR, DriftSpacerUpperSquareCutout_oR, DriftSpacerUpperSquareCutout_thick/2, opendeg, closedeg);
	G4Tubs *DriftSpacerLowerSquareCutout= new G4Tubs("DriftSpacerLowerSquareCutout",	DriftSpacerLowerSquareCutout_iR, DriftSpacerLowerSquareCutout_oR, DriftSpacerLowerSquareCutout_thick/2, opendeg, closedeg);

	G4ThreeVector DriftSpacerUpperSquareCutout_V	(0.0, 0.0, 0.0);
	G4ThreeVector DriftSpacerLowerSquareCutout_V	(0.0, 0.0, 0.0 );
	
	G4Transform3D DriftSpacerUpperSquareCutout_T	(ZeroRot, DriftSpacerUpperSquareCutout_V);
	G4Transform3D DriftSpacerLowerSquareCutout_T	(ZeroRot, DriftSpacerLowerSquareCutout_V);

	G4SubtractionSolid *DriftSpacerUpperSquare1 = new G4SubtractionSolid("DriftSpacerUpperSquare1",	DriftSpacerUpperSquare,	DriftSpacerUpperSquareCutout,	DriftSpacerUpperSquareCutout_T);
	G4SubtractionSolid *DriftSpacerLowerSquare1 = new G4SubtractionSolid("DriftSpacerLowerSquare1",	DriftSpacerLowerSquare,	DriftSpacerLowerSquareCutout,	DriftSpacerLowerSquareCutout_T);
	
	// Bottom PMT Holder
	G4double BottomHolderTube_iR				= 0.0 *mm;
	G4double BottomHolderTube_oR				= 60.0 *mm;
	G4double BottomHolderTube_height			= 33.0 *mm;
	
	G4double BottomHolderUpperTube_iR			= 17.5 *mm;
	G4double BottomHolderUpperTube_oR			= 60.0 *mm;
	G4double BottomHolderUpperTube_height		= 5.0 *mm;

	G4double BottomHolderMiddleTube_iR			= 28.5 *mm;
	G4double BottomHolderMiddleTube_oR			= 60.0 *mm;
	G4double BottomHolderMiddleTube_height		= 26.0 *mm;

	G4double BottomHolderLowerTube_iR			= 28.5 *mm;
	G4double BottomHolderLowerTube_oR			= 60.0 *mm;
	G4double BottomHolderLowerTube_height		= 2.0 *mm;
	
	G4Tubs *BottomHolderUpperTube 	= new G4Tubs("BottomHolderUpperTube", BottomHolderUpperTube_iR, BottomHolderUpperTube_oR, BottomHolderUpperTube_height/2, opendeg, closedeg);
	G4Tubs *BottomHolderMiddleTube	= new G4Tubs("BottomHolderMiddleTube", BottomHolderMiddleTube_iR, BottomHolderMiddleTube_oR, BottomHolderMiddleTube_height/2, opendeg, closedeg);
	G4Tubs *BottomHolderLowerTube	= new G4Tubs("BottomHolderLowerTube", BottomHolderLowerTube_iR, BottomHolderLowerTube_oR, BottomHolderLowerTube_height/2, opendeg, closedeg);
	
	// Bottom PMT Clamp
	G4double BottomClamp_height				= 9.0 *mm;
	G4double BottomClamp_iR					= 22.5 *mm;

	G4double BottomClampUpperPart_height	= 2.0 *mm;
	G4double BottomClampUpperPart_oR		= 30.0 *mm;
	
	G4double BottomClampLowerPart_height	= 2.0 *mm;
	G4double BottomClampLowerPart_oR		= 23.5 *mm;
	
	G4double BottomClampMiddlePart_height	= 5.0 *mm;
	G4double BottomClampMiddlePart_oR		= 23.5 *mm;
	
	G4Tubs *BottomClampMiddlePart	= new G4Tubs("BottomClampMiddlePart", BottomClamp_iR, BottomClampMiddlePart_oR, BottomClampMiddlePart_height/2, opendeg, closedeg);
	G4Tubs *BottomClampUpperPart	= new G4Tubs("BottomClampUpperPart", BottomClamp_iR, BottomClampUpperPart_oR, BottomClampUpperPart_height/2, opendeg, closedeg);
	G4Tubs *BottomClampLowerPart	= new G4Tubs("BottomClampLowerPart", BottomClamp_iR, BottomClampLowerPart_oR, BottomClampLowerPart_height/2, opendeg, closedeg);
	
	//--------------	Filler, with cut filled with GXe and LXe
	G4double Filler_oR				= 73.75 *mm;
	G4double Filler_height			= 65.0 *mm;	
	G4double FillerUpperPart_iR		= 60.5 *mm;
	G4double FillerUpperPart_oR		= Filler_oR;
	G4double FillerUpperPart_height	= 63.0 *mm;
	G4double FillerLowerPart_iR		= 58.0 *mm;
	G4double FillerLowerPart_oR		= Filler_oR;
	G4double FillerLowerPart_height	= 2.0 *mm;
	G4double FillerCut_iR			= 63.125 *mm;
	G4double FillerCut_oR			= 71.125 *mm;
	G4double FillerCut_depth		= 50.0 *mm;

	G4Tubs *FillerUpperPart	= new G4Tubs("FillerUpperPart", FillerUpperPart_iR, FillerUpperPart_oR, FillerUpperPart_height/2, opendeg, closedeg);
	G4Tubs *FillerLowerPart	= new G4Tubs("FillerLowerPart", FillerLowerPart_iR, FillerLowerPart_oR, FillerLowerPart_height/2, opendeg, closedeg);

	G4ThreeVector FillerLowerPart_V		(0.0, 0.0, -FillerUpperPart_height/2-FillerLowerPart_height/2 );	
	G4Transform3D FillerLowerPart_T		(ZeroRot, FillerLowerPart_V);

	G4UnionSolid *Filler1	= new G4UnionSolid(	"Filler1",	FillerUpperPart,	FillerLowerPart,	FillerLowerPart_T);

	G4Tubs *FillerCut	= new G4Tubs("FillerCut", FillerCut_iR, FillerCut_oR, FillerCut_depth/2, opendeg, closedeg);

	G4ThreeVector FillerCut_V	(0.0, 0.0,	Filler_height/2 - FillerCut_depth/2);
	G4Transform3D FillerCut_T	(ZeroRot, FillerCut_V);

	G4SubtractionSolid *Filler = new G4SubtractionSolid( "Filler",	Filler1, FillerCut,	FillerCut_T);

	// Grids
	G4double GridFrame_thick		= 3.0 *mm;
	G4double GridFrame_side			= 84.0 *mm;
	G4double GridFrameCutout_thick	= GridFrame_thick;
	G4double GridFrameCutout_side	= 59.0 *mm;

	G4Box *GridFrameBox		= new G4Box("GridFrameBox", GridFrame_side/2, GridFrame_side/2, GridFrame_thick/2);
	G4Box *GridFrameCutout	= new G4Box("GridFrameCutout", GridFrameCutout_side/2, GridFrameCutout_side/2, GridFrameCutout_thick/2);

	G4ThreeVector GridFrameCutout_V	(0.0, 0.0,	0.0 );	
	G4Transform3D GridFrameCutout_T	(ZeroRot, GridFrameCutout_V);

	G4SubtractionSolid *GridFrame	= new G4SubtractionSolid( "GridFrame",	GridFrameBox,	GridFrameCutout,	GridFrameCutout_T);
	// distance between the grid wires = 12/3, 29 wires, thick 3
	
	
	// Photomultiplier Tube
	G4double PMTcasing_height	= 32.5 *mm;
	G4double PMTcasing_iR		= 0.0 *mm;
	G4double PMTcasing_oR		= 28.5 *mm;
	G4double PMTcasing_thick	= 0.5 *mm;
	G4double PMTinterior_height	= PMTcasing_height-PMTcasing_thick;
	G4double PMTinterior_iR		= 0.0 *mm;
	G4double PMTinterior_oR		= PMTcasing_oR-PMTcasing_thick;
	G4double PMTcathode_thick	= 0.5 *mm;
	G4double PMTcathode_iR		= 0.0 *mm;
	G4double PMTcathode_oR		= 22.5 *mm;
	G4double PMTwindow_thick	= 1.5 *mm;
	G4double PMTwindow_iR		= 0.0 *mm;
	G4double PMTwindow_oR		= 25.25 *mm;

	G4Tubs *PMTcasing	= new G4Tubs("PMTcasing",	PMTcasing_iR, PMTcasing_oR, PMTcasing_height/2, opendeg, closedeg);
	G4Tubs *PMTinterior	= new G4Tubs("PMTinterior",	PMTinterior_iR, PMTinterior_oR, PMTinterior_height/2, opendeg, closedeg);
	G4Tubs *PMTcathode 	= new G4Tubs("PMTcathode",	PMTcathode_iR, PMTcathode_oR, PMTcathode_thick/2, opendeg, closedeg);
	G4Tubs *PMTwindow	= new G4Tubs("PMTwindow",	PMTwindow_iR, PMTwindow_oR, PMTwindow_thick/2, opendeg, closedeg);

	// PMT base
	G4double PMTbase_thick	= 5.0 *mm;
	G4double PMTbase_iR		= 7.5 *mm;
	G4double PMTbase_oR		= 20.0 *mm;

	G4Tubs *PMTbase	= new G4Tubs("PMTbase", PMTbase_iR, PMTbase_oR, PMTbase_thick/2, opendeg, closedeg);	


	////////////////////////////////		
	// NaI detector
	G4double NaI_height 		= 222.25 *mm;
	
	G4double NaI_Crystal_iD 	= 0.0 *mm;
	G4double NaI_Crystal_oD 	= 76.2 *mm;
	G4double NaI_Crystal_height = 76.2 *mm;
	
	G4double NaI_CrystalHousing_thick		= 0.508 *mm;
	G4double NaI_CrystalHousingTop_oD		= 80.72 *mm;
	G4double NaI_CrystalHousingTop_iD		= NaI_CrystalHousingTop_oD - 2*NaI_CrystalHousing_thick;
	G4double NaI_CrystalHousing_height		= NaI_Crystal_height + NaI_CrystalHousing_thick;
	G4double NaI_CrystalHousingTop_height	= NaI_Crystal_height;
	G4double NaI_CrystalHousingBot_height	= NaI_CrystalHousing_thick;
	G4double NaI_CrystalHousingBot_iD		= 0.0 *mm;
	G4double NaI_CrystalHousingBot_oD		= NaI_CrystalHousingTop_oD;

	G4double NaI_LightShield_thick			= 0.635 *mm;
	G4double NaI_LightShield_height			= NaI_height - NaI_CrystalHousing_height;
	G4double NaI_LightShieldBot_oD			= 82.55 *mm;
	G4double NaI_LightShieldBot_iD			= NaI_LightShieldBot_oD - 2*NaI_LightShield_thick;
	G4double NaI_LightShieldBot_height		= 30.00 *mm;

	G4double NaI_LightShieldTop_oD			= 58.72 *mm;
	G4double NaI_LightShieldTop_iD			= NaI_LightShieldTop_oD - 2*NaI_LightShield_thick;

	G4double NaI_LightShieldMid_oRtop		= NaI_LightShieldTop_oD/2;
	G4double NaI_LightShieldMid_iRtop		= NaI_LightShieldMid_oRtop - NaI_LightShield_thick;
	G4double NaI_LightShieldMid_oRbot		= NaI_LightShieldBot_oD/2;
	G4double NaI_LightShieldMid_iRbot		= NaI_LightShieldMid_oRbot - NaI_LightShield_thick;
	G4double NaI_LightShieldMid_height		= 15.0 *mm;

	G4double NaI_LightShieldTop_height		= NaI_LightShield_height - NaI_LightShieldBot_height - NaI_LightShieldMid_height;
	// vacuum inside the light shield
/*	G4double NaI_LightShieldVac_height		= NaI_height - NaI_CrystalHousing_height;
	G4double NaI_LightShieldVacBot_oD		= NaI_LightShieldBot_iD;
	G4double NaI_LightShieldVacBot_iD		= 0.0 *mm;
	G4double NaI_LightShieldVacBot_height	= 30.00 *mm;

	G4double NaI_LightShieldVacTop_oD		= NaI_LightShieldTop_iD;
	G4double NaI_LightShieldVacTop_iD		= 0.0 *mm;

	G4double NaI_LightShieldVacMid_oRtop	= NaI_LightShieldMid_iRtop;
	G4double NaI_LightShieldVacMid_iRtop	= 0.0 *mm;
	G4double NaI_LightShieldVacMid_oRbot	= NaI_LightShieldMid_iRbot;
	G4double NaI_LightShieldVacMid_iRbot	= 0.0 *mm;
	G4double NaI_LightShieldVacMid_height	= 19.4 *mm;

	G4double NaI_LightShieldVacTop_height	= NaI_LightShield_height - NaI_LightShieldBot_height - NaI_LightShieldMid_height;
*/
	// NaI PMT	
	G4double NaI_PMTcasing_thick		= 1.0 *mm;
	G4double NaI_PMTcasing_height		= 123.0 *mm;
	G4double NaI_PMTcasingBot_height	= 30.0 *mm;
	G4double NaI_PMTcasingMid_height	= 19.4 *mm;
	G4double NaI_PMTcasingTop_height	= NaI_PMTcasing_height - NaI_PMTcasingBot_height - NaI_PMTcasingMid_height;
	G4double NaI_PMTcasingBot_oD		= 78.0 *mm;
	G4double NaI_PMTcasingTop_oD		= 51.5 *mm;
	G4double NaI_PMTcasingMid_oDbot		= NaI_PMTcasingBot_oD;
	G4double NaI_PMTcasingMid_oDtop		= NaI_PMTcasingTop_oD;
	G4double NaI_PMTcasingBot_iD		= NaI_PMTcasingBot_oD - 2*NaI_PMTcasing_thick;
	G4double NaI_PMTcasingTop_iD		= NaI_PMTcasingTop_oD - 2*NaI_PMTcasing_thick;
	G4double NaI_PMTcasingMid_iDbot		= NaI_PMTcasingBot_iD;
	G4double NaI_PMTcasingMid_iDtop		= NaI_PMTcasingTop_iD;
		
	G4Tubs *NaI_Crystal				= new G4Tubs("NaI_Crystal",	NaI_Crystal_iD/2, NaI_Crystal_oD/2, NaI_Crystal_height/2, opendeg, closedeg);
	G4Tubs *NaI_CrystalHousingTop	= new G4Tubs("NaI_CrystalHousingTop",	NaI_CrystalHousingTop_iD/2, NaI_CrystalHousingTop_oD/2, NaI_CrystalHousingTop_height/2, opendeg, closedeg);
	G4Tubs *NaI_CrystalHousingBot	= new G4Tubs("NaI_CrystalHousingBot",	NaI_CrystalHousingBot_iD/2, NaI_CrystalHousingBot_oD/2, NaI_CrystalHousingBot_height/2, opendeg, closedeg);

	G4Tubs *NaI_LightShieldBot		= new G4Tubs("NaI_LightShieldBot",	NaI_LightShieldBot_iD/2, NaI_LightShieldBot_oD/2, NaI_LightShieldBot_height/2, opendeg, closedeg);
	G4Tubs *NaI_LightShieldTop		= new G4Tubs("NaI_LightShieldTop",	NaI_LightShieldTop_iD/2, NaI_LightShieldTop_oD/2, NaI_LightShieldTop_height/2, opendeg, closedeg);
	G4Cons *NaI_LightShieldMid		= new G4Cons("NaI_LightShieldMid",	NaI_LightShieldMid_iRbot, NaI_LightShieldMid_oRbot, NaI_LightShieldMid_iRtop, NaI_LightShieldMid_oRtop, NaI_LightShieldMid_height/2, opendeg, closedeg);
	//G4ThreeVector NaI_LightShieldMid_V	(0.0, 0.0, NaI_LightShieldBot_height/2+NaI_LightShieldMid_height/2);
	//G4Transform3D NaI_LightShieldMid_T	(ZeroRot, NaI_LightShieldMid_V);
	//G4UnionSolid *NaI_LightShield1	= new G4UnionSolid( "NaI_LightShield1",	NaI_LightShieldBot, NaI_LightShieldMid,	NaI_LightShieldMid_T);
	//G4ThreeVector NaI_LightShieldTop_V	(0.0, 0.0, (NaI_LightShieldBot_height+NaI_LightShieldMid_height)/2+NaI_LightShieldTop_height/2);
	//G4Transform3D NaI_LightShieldTop_T	(ZeroRot, NaI_LightShieldTop_V);
	//G4UnionSolid *NaI_LightShield	= new G4UnionSolid( "NaI_LightShield",	NaI_LightShield1,	NaI_LightShieldTop,	NaI_LightShieldTop_T);

/*	G4Tubs *NaI_LightShieldVacBot	= new G4Tubs("NaI_LightShieldVacBot",	NaI_LightShieldVacBot_iD/2, NaI_LightShieldVacBot_oD/2, NaI_LightShieldVacBot_height/2, opendeg, closedeg);
	G4Tubs *NaI_LightShieldVacTop	= new G4Tubs("NaI_LightShieldVacTop",	NaI_LightShieldVacTop_iD/2, NaI_LightShieldVacTop_oD/2, NaI_LightShieldVacTop_height/2, opendeg, closedeg);
	G4Cons *NaI_LightShieldVacMid	= new G4Cons("NaI_LightShieldVacMid",	NaI_LightShieldVacMid_iRbot, NaI_LightShieldVacMid_oRbot, NaI_LightShieldVacMid_iRtop, NaI_LightShieldVacMid_oRtop, NaI_LightShieldVacMid_height/2, opendeg, closedeg);
	G4ThreeVector NaI_LightShieldVacMid_V	(0.0, 0.0, NaI_LightShieldVacBot_height/2+NaI_LightShieldVacMid_height/2);
	G4Transform3D NaI_LightShieldVacMid_T	(ZeroRot, NaI_LightShieldVacMid_V);
	G4UnionSolid *NaI_LightShieldVac1	= new G4UnionSolid( "NaI_LightShieldVac1", NaI_LightShieldVacBot, NaI_LightShieldVacMid, NaI_LightShieldVacMid_T);
	G4ThreeVector NaI_LightShieldVacTop_V	(0.0, 0.0, (NaI_LightShieldVacBot_height+NaI_LightShieldVacMid_height)/2+NaI_LightShieldVacTop_height/2);
	G4Transform3D NaI_LightShieldVacTop_T	(ZeroRot, NaI_LightShieldVacTop_V);
	G4UnionSolid *NaI_LightShieldVac	= new G4UnionSolid( "NaI_LightShieldVac", NaI_LightShieldVac1, NaI_LightShieldVacTop, NaI_LightShieldVacTop_T);
*/
/*	G4Tubs *NaI_PMTcasingBot	= new G4Tubs("NaI_PMTcasingBot", NaI_PMTcasingBot_iD/2, NaI_PMTcasingBot_oD/2, NaI_PMTcasingBot_height/2, opendeg, closedeg);
	G4Tubs *NaI_PMTcasingTop	= new G4Tubs("NaI_PMTcasingTop", NaI_PMTcasingTop_iD/2, NaI_PMTcasingTop_oD/2, NaI_PMTcasingTop_height/2, opendeg, closedeg);
	G4Cons *NaI_PMTcasingMid	= new G4Cons("NaI_PMTcasingMid", NaI_PMTcasingMid_iDbot/2, NaI_PMTcasingMid_oDbot/2, NaI_PMTcasingMid_iDtop/2, NaI_PMTcasingMid_oDtop/2, NaI_PMTcasingMid_height/2, opendeg, closedeg);
	G4ThreeVector NaI_PMTcasingMid_V	(0.0, 0.0, NaI_PMTcasingBot_height/2+NaI_PMTcasingMid_height/2);
	G4Transform3D NaI_PMTcasingMid_T	(ZeroRot, NaI_PMTcasingMid_V);
	G4UnionSolid *NaI_PMTcasing1	= new G4UnionSolid( "NaI_PMTcasing1",	NaI_PMTcasingBot, NaI_PMTcasingMid,	NaI_PMTcasingMid_T);
	G4ThreeVector NaI_PMTcasingTop_V	(0.0, 0.0, (NaI_PMTcasingBot_height+NaI_PMTcasingMid_height)/2+NaI_PMTcasingTop_height/2);
	G4Transform3D NaI_PMTcasingTop_T	(ZeroRot, NaI_PMTcasingTop_V);
	G4UnionSolid *NaI_PMTcasing	= new G4UnionSolid( "NaI_PMTcasing",	NaI_PMTcasing1,	NaI_PMTcasingTop,	NaI_PMTcasingTop_T);
*/
	////////////////////////////////////////////////////////////////////////////
	// STAND FOR NaI detector	
	G4double NaIstand_profile_thick 		= 4.5 *mm;	
	G4double NaIstand_profileOut_out 		= 40.0 *mm;
	G4double NaIstand_profileOut_in 		= NaIstand_profileOut_out - 2*NaIstand_profile_thick;
	G4double NaIstand_profileIn_in	 		= 13.7 *mm;
	G4double NaIstand_profileIn_out 		= NaIstand_profileIn_in + 2*NaIstand_profile_thick;
	G4double NaIstand_shortProfile_length 	= 320.0 *mm;
	G4double NaIstand_longProfile_height 	= 766.0 *mm;

	G4double NaIstand_Plate_thick 			= 10.0 *mm;
	G4double NaIstand_Plate_side 			= 400.0 *mm;
	
	G4double NaIstand_CradleBar_height 		= 670.0 *mm;
	G4double NaIstand_CradleBar_iD 			= 0.0 *mm;
	G4double NaIstand_CradleBar_oD 			= 20.0 *mm;

	G4double NaIstand_CradlePanel_thick 	= 10.0 *mm;
	G4double NaIstand_CradlePanel_length 	= 380.0 *mm;
	G4double NaIstand_CradlePanel_width 	= 20.0 *mm;
		
	G4Box *NaIstand_VertBarOut_out 		= new G4Box("NaIstand_VertBarOut_out", NaIstand_profileOut_out/2, NaIstand_profileOut_out/2, NaIstand_longProfile_height/2);
	G4Box *NaIstand_VertBarOut_in 		= new G4Box("NaIstand_VertBarOut_in", NaIstand_profileOut_in/2, NaIstand_profileOut_in/2, NaIstand_longProfile_height/2);
	G4ThreeVector NaIstand_VertBarOut_V (0.0, 0.0,	0.0);
	G4Transform3D NaIstand_VertBarOut_T (ZeroRot, NaIstand_VertBarOut_V);	
	G4SubtractionSolid *NaIstand_VertBarOut	= new G4SubtractionSolid( "NaIstand_VertBarOut", NaIstand_VertBarOut_out, NaIstand_VertBarOut_in, NaIstand_VertBarOut_T);

	G4Box *NaIstand_VertBarIn_out 		= new G4Box("NaIstand_VertBarIn_out", NaIstand_profileIn_out/2, NaIstand_profileIn_out/2, NaIstand_longProfile_height/2);
	G4Box *NaIstand_VertBarIn_in 		= new G4Box("NaIstand_VertBarIn_in", NaIstand_profileIn_in/2, NaIstand_profileIn_in/2, NaIstand_longProfile_height/2);
	G4ThreeVector NaIstand_VertBarIn_V (0.0, 0.0,	0.0);
	G4Transform3D NaIstand_VertBarIn_T (ZeroRot, NaIstand_VertBarIn_V);	
	G4SubtractionSolid *NaIstand_VertBarIn	= new G4SubtractionSolid( "NaIstand_VertBarIn", NaIstand_VertBarIn_out, NaIstand_VertBarIn_in, NaIstand_VertBarIn_T);

	G4Box *NaIstand_HorzBarXout_out 	= new G4Box("NaIstand_HorzBarXout_out", NaIstand_shortProfile_length/2, NaIstand_profileOut_out/2, NaIstand_profileOut_out/2);
	G4Box *NaIstand_HorzBarXout_in 		= new G4Box("NaIstand_HorzBarXout_in", NaIstand_shortProfile_length/2, NaIstand_profileOut_in/2, NaIstand_profileOut_in/2);
	G4ThreeVector NaIstand_HorzBarXout_V (0.0, 0.0,	0.0);
	G4Transform3D NaIstand_HorzBarXout_T (ZeroRot, NaIstand_HorzBarXout_V);	
	G4SubtractionSolid *NaIstand_HorzBarXout	= new G4SubtractionSolid( "NaIstand_HorzBarXout", NaIstand_HorzBarXout_out, NaIstand_HorzBarXout_in, NaIstand_HorzBarXout_T);

	G4Box *NaIstand_HorzBarXin_out 		= new G4Box("NaIstand_HorzBarXin_out", NaIstand_shortProfile_length/2, NaIstand_profileOut_out/2, NaIstand_profileOut_out/2);
	G4Box *NaIstand_HorzBarXin_in 		= new G4Box("NaIstand_HorzBarXin_in", NaIstand_shortProfile_length/2, NaIstand_profileOut_in/2, NaIstand_profileIn_in/2);
	G4ThreeVector NaIstand_HorzBarXin_V (0.0, 0.0,	0.0);
	G4Transform3D NaIstand_HorzBarXin_T (ZeroRot, NaIstand_HorzBarXin_V);	
	G4SubtractionSolid *NaIstand_HorzBarXin	= new G4SubtractionSolid( "NaIstand_HorzBarXin", NaIstand_HorzBarXin_out, NaIstand_HorzBarXin_in, NaIstand_HorzBarXin_T);

	G4Box *NaIstand_HorzBarYout_out 	= new G4Box("NaIstand_HorzBarYout_out", NaIstand_profileOut_out/2, NaIstand_shortProfile_length/2, NaIstand_profileOut_out/2);
	G4Box *NaIstand_HorzBarYout_in 		= new G4Box("NaIstand_HorzBarYout_in", NaIstand_profileOut_in/2, NaIstand_shortProfile_length/2, NaIstand_profileOut_in/2);
	G4ThreeVector NaIstand_HorzBarYout_V (0.0, 0.0,	0.0);
	G4Transform3D NaIstand_HorzBarYout_T (ZeroRot, NaIstand_HorzBarYout_V);	
	G4SubtractionSolid *NaIstand_HorzBarYout	= new G4SubtractionSolid( "NaIstand_HorzBarYout", NaIstand_HorzBarYout_out, NaIstand_HorzBarYout_in, NaIstand_HorzBarYout_T);

	G4Box *NaIstand_HorzBarYin_out 		= new G4Box("NaIstand_HorzBarYin_out", NaIstand_profileOut_out/2, NaIstand_shortProfile_length/2, NaIstand_profileOut_out/2);
	G4Box *NaIstand_HorzBarYin_in 		= new G4Box("NaIstand_HorzBarYin_in", NaIstand_profileOut_in/2, NaIstand_shortProfile_length/2, NaIstand_profileIn_in/2);
	G4ThreeVector NaIstand_HorzBarYin_V (0.0, 0.0,	0.0);
	G4Transform3D NaIstand_HorzBarYin_T (ZeroRot, NaIstand_HorzBarYin_V);	
	G4SubtractionSolid *NaIstand_HorzBarYin	= new G4SubtractionSolid( "NaIstand_HorzBarYin", NaIstand_HorzBarYin_out, NaIstand_HorzBarYin_in, NaIstand_HorzBarYin_T);

	G4Box *NaIstand_Plate 				= new G4Box("NaIstand_Plate", NaIstand_Plate_side/2, NaIstand_Plate_side/2, NaIstand_Plate_thick/2);

	G4Tubs *NaIstand_CradleBar			= new G4Tubs("NaIstand_CradleBar",	NaIstand_CradleBar_iD/2, NaIstand_CradleBar_oD/2, NaIstand_CradleBar_height/2, opendeg, closedeg);
	G4Box  *NaIstand_CradlePanel 		= new G4Box("NaIstand_CradlePanel", NaIstand_CradlePanel_length/2, NaIstand_CradlePanel_width/2, NaIstand_CradlePanel_thick/2);

	////////////////////////////////////////////////////////////////////////////
	// LEAD CASTLE FOR NaI detector	
	G4double NaIcastle_thick 		= 5.0 *cm;	
	G4double NaIcastle_width 		= 20.0 *cm;
	G4double NaIcastle_height 		= 20.0 *cm;
	G4double NaIcastle_length	 	= 45.0 *cm;

	G4double NaIcastleCavity_width 	= NaIcastle_width - 2*NaIcastle_thick;
	G4double NaIcastleCavity_height = NaIcastle_height - 2*NaIcastle_thick;
	G4double NaIcastleCavity_length = NaIcastle_length - NaIcastle_thick;

	G4Box  *NaIcastleBox 	= new G4Box("NaIcastleBox", NaIcastle_width/2, NaIcastle_length/2, NaIcastle_height/2);
	G4Box  *NaIcastleCavity = new G4Box("NaIcastleCavity", NaIcastleCavity_width/2, NaIcastleCavity_length/2, NaIcastleCavity_height/2);

	G4ThreeVector NaIcastleCavity_V (0.0, NaIcastle_length/2-NaIcastleCavity_length/2, 0.0);
	G4Transform3D NaIcastleCavity_T (ZeroRot, NaIcastleCavity_V);
	G4SubtractionSolid *NaIcastle	= new G4SubtractionSolid( "NaIcastle", NaIcastleBox, NaIcastleCavity, NaIcastleCavity_T);

	////////////////////////////////////////////////////////////////////////////
	// Collimator	
	G4double CollimatorBox_oD 		= 132.0 *mm;	
	G4double CollimatorBox_iD 		= 0.0;
	G4double CollimatorBox_height 	= 100.0 *mm;
	G4double CollimatorHole_oD 		= 6.0 *mm;	
	G4double CollimatorHole_iD 		= 0.0;
	G4double CollimatorHole_height 	= 50.0 *mm;
	G4Tubs *CollimatorBox	= new G4Tubs("ColllimatorBox",	CollimatorBox_iD/2, CollimatorBox_oD/2, CollimatorBox_height/2, opendeg, closedeg);
	G4Tubs *CollimatorHole	= new G4Tubs("ColllimatorHole",	CollimatorHole_iD/2, CollimatorHole_oD/2, CollimatorHole_height/2, opendeg, closedeg);
	G4ThreeVector CollimatorHole_V (0.0, 0.0, CollimatorBox_height/2-CollimatorHole_height/2);
	G4Transform3D CollimatorHole_T (ZeroRot, CollimatorHole_V);
	G4SubtractionSolid *Collimator	= new G4SubtractionSolid( "Collimator", CollimatorBox, CollimatorHole, CollimatorHole_T);

	////////////////////////////////////////////////////////////////////////////
	// LEAD CHANNEL
	G4double LeadBrick_height 	= 10. *cm;
	G4double LeadBrick_length 	= 20. *cm;
	G4double LeadBrick_thick 	= 5. *cm;
	
	G4double LeadChannel_thick 	= 5.0 *cm;	
	G4double LeadChannel_width 	= 15.0 *cm;
	G4double LeadChannel_height 	= 20.0 *cm;
	G4double LeadChannel_length	= 80.0 *cm;

	G4double LeadChannelCavitySide_width 	= 5. *cm;
	G4double LeadChannelCavitySide_height 	= 10. *cm;
	G4double LeadChannelCavitySide_length 	= 20. *cm;

	G4double LeadChannelCavityCenter_width 	= 3*cm;
	G4double LeadChannelCavityCenter_height = 10. *cm;
	G4double LeadChannelCavityCenter_length = 40. *cm;

	G4Box  *LeadChannelBox 			= new G4Box("LeadChannelBox", LeadChannel_width/2, LeadChannel_length/2, LeadChannel_height/2);
	G4Box  *LeadChannelCavitySide 	= new G4Box("LeadChannelCavitySide", LeadChannelCavitySide_width/2, LeadChannelCavitySide_length/2, LeadChannelCavitySide_height/2);
	G4Box  *LeadChannelCavityCenter = new G4Box("LeadChannelCavityCenter", LeadChannelCavityCenter_width/2, LeadChannelCavityCenter_length/2, LeadChannelCavityCenter_height/2);

	G4ThreeVector LeadChannelCavitySide1_V (0.0, LeadChannel_length/2-LeadChannelCavitySide_length/2, 0.0);
	G4Transform3D LeadChannelCavitySide1_T (ZeroRot, LeadChannelCavitySide1_V);
	G4ThreeVector LeadChannelCavitySide2_V (0.0, -LeadChannel_length/2+LeadChannelCavitySide_length/2, 0.0);
	G4Transform3D LeadChannelCavitySide2_T (ZeroRot, LeadChannelCavitySide2_V);
	G4ThreeVector LeadChannelCavityCenter_V (0.0, 0.0, 0.0);
	G4Transform3D LeadChannelCavityCenter_T (ZeroRot, LeadChannelCavityCenter_V);

	G4SubtractionSolid *LeadChannel2	= new G4SubtractionSolid( "LeadChannel2", LeadChannelBox, LeadChannelCavitySide, LeadChannelCavitySide1_T);
	G4SubtractionSolid *LeadChannel1	= new G4SubtractionSolid( "LeadChannel1", LeadChannel2, LeadChannelCavitySide, LeadChannelCavitySide2_T);
	G4SubtractionSolid *LeadChannel		= new G4SubtractionSolid( "LeadChannel", LeadChannel1, LeadChannelCavityCenter, LeadChannelCavityCenter_T);

	////////////////////////////////////////////////////////////////////////////
	// LEAD APERTURE
	G4double Aperture_thick 	= 5.0 *cm;
	G4double Aperture_width 	= 10. *cm;	//9.5 *cm;
	G4double Aperture_height 	= 10. *cm;	//9.5 *cm;

	G4double ApertureHole_thick = Aperture_thick;
	G4double ApertureHole_iD 	= 0.0;
	G4double ApertureHole_oD 	= 3.0 *cm;

	G4Box  *ApertureBox 	= new G4Box("ApertureBox", Aperture_width/2, Aperture_thick/2, Aperture_height/2);
	G4Tubs  *ApertureHole 	= new G4Tubs("ApertureHole", ApertureHole_iD/2, ApertureHole_oD/2, ApertureHole_thick/2, opendeg, closedeg);

	G4ThreeVector ApertureHole_V (0.0, 0.0, 0.0);
	G4Transform3D ApertureHole_T (RotationXPlus90, ApertureHole_V);
	G4SubtractionSolid *Aperture	= new G4SubtractionSolid( "Aperture", ApertureBox, ApertureHole, ApertureHole_T);

	G4double ApertureFrame_thick 	= 5.0 *cm;
	G4double ApertureFrame_width 	= 16.3 *cm;
	G4double ApertureFrame_height 	= 16.3 *cm;
	G4double ApertureFrame_profile 	= 15. *mm;
	
	G4double ApertureFrameHole_thick 	= Aperture_thick;
	G4double ApertureFrameHole_width 	= ApertureFrame_width - 2*ApertureFrame_profile;
	G4double ApertureFrameHole_height 	= ApertureFrame_height - 2*ApertureFrame_profile;

	G4Box  *ApertureFrameBox 	= new G4Box("ApertureFrameBox", ApertureFrame_width/2, ApertureFrame_thick/2, ApertureFrame_height/2);
	G4Box  *ApertureFrameHole 	= new G4Box("ApertureFrameHole", ApertureFrameHole_width/2, ApertureFrameHole_thick/2, ApertureFrameHole_height/2);

	G4ThreeVector ApertureFrameHole_V (0.0, 0.0, 0.0);
	G4Transform3D ApertureFrameHole_T (ZeroRot, ApertureFrameHole_V);
	G4SubtractionSolid *ApertureFrame	= new G4SubtractionSolid( "ApertureFrame", ApertureFrameBox, ApertureFrameHole, ApertureFrameHole_T);
	
	//==============================================================================================
	//==============	Logical Volumes (declared in 'XuerichDetectorGeometry.hh')	================
	//==============================================================================================
	G4LogicalVolume *Laboratory_log					= new G4LogicalVolume( Laboratory, 				Vacuum, 	"Laboratory_log");
	// stand for xuerich
	G4LogicalVolume	*XuStand_VertBar_LeftFront_log	= new G4LogicalVolume( XuStand_VertBar,			metalAl,	"XuStand_VertBar_LeftFront_log");
	G4LogicalVolume	*XuStand_VertBar_LeftBack_log	= new G4LogicalVolume( XuStand_VertBar,			metalAl,	"XuStand_VertBar_LeftBack_log");
	G4LogicalVolume	*XuStand_VertBar_RightFront_log	= new G4LogicalVolume( XuStand_VertBar,			metalAl,	"XuStand_VertBar_RightFront_log");
	G4LogicalVolume	*XuStand_VertBar_RightBack_log	= new G4LogicalVolume( XuStand_VertBar,			metalAl,	"XuStand_VertBar_RightBack_log");

	G4LogicalVolume	*XuStand_HorzBar_BotLeft_log	= new G4LogicalVolume( XuStand_HorzBarShortY,	metalAl,	"XuStand_HorzBar_BotLeft_log");
	G4LogicalVolume	*XuStand_HorzBar_BotRight_log	= new G4LogicalVolume( XuStand_HorzBarShortY,	metalAl,	"XuStand_HorzBar_BotRight_log");
	G4LogicalVolume	*XuStand_HorzBar_BotFront_log	= new G4LogicalVolume( XuStand_HorzBarLongX,	metalAl,	"XuStand_HorzBar_BotFront_log");
	G4LogicalVolume	*XuStand_HorzBar_BotBack_log	= new G4LogicalVolume( XuStand_HorzBarLongX,	metalAl,	"XuStand_HorzBar_BotBack_log");

	G4LogicalVolume	*XuStand_HorzBar_MidLeft_log	= new G4LogicalVolume( XuStand_HorzBarShortY,	metalAl,	"XuStand_HorzBar_MidLeft_log");
	G4LogicalVolume	*XuStand_HorzBar_MidRight_log	= new G4LogicalVolume( XuStand_HorzBarShortY,	metalAl,	"XuStand_HorzBar_MidRight_log");
	G4LogicalVolume	*XuStand_HorzBar_MidFront_log	= new G4LogicalVolume( XuStand_HorzBarShortX,	metalAl,	"XuStand_HorzBar_MidFront_log");
	G4LogicalVolume	*XuStand_HorzBar_MidBack_log	= new G4LogicalVolume( XuStand_HorzBarShortX,	metalAl,	"XuStand_HorzBar_MidBack_log");

	G4LogicalVolume	*XuStand_HorzBar_TopLeft_log	= new G4LogicalVolume( XuStand_HorzBarLongY,	metalAl,	"XuStand_HorzBar_TopLeft_log");
	G4LogicalVolume	*XuStand_HorzBar_TopRight_log	= new G4LogicalVolume( XuStand_HorzBarLongY,	metalAl,	"XuStand_HorzBar_TopRight_log");
	G4LogicalVolume	*XuStand_HorzBar_TopFront_log	= new G4LogicalVolume( XuStand_HorzBarShortX,	metalAl,	"XuStand_HorzBar_TopFront_log");
	G4LogicalVolume	*XuStand_HorzBar_TopBack_log	= new G4LogicalVolume( XuStand_HorzBarShortX,	metalAl,	"XuStand_HorzBar_TopBack_log");

	G4LogicalVolume	*XuStand_BottomPlate_log		= new G4LogicalVolume( XuStand_BottomPlate,		metalAl,	"XuStand_BottomPlate_log");

	G4LogicalVolume	*XuStand_WingFront_log			= new G4LogicalVolume( XuStand_Wing,			metalAl,	"XuStand_WingFront_log");
	G4LogicalVolume	*XuStand_WingBack_log			= new G4LogicalVolume( XuStand_Wing,			metalAl,	"XuStand_WingBack_log");

	G4LogicalVolume	*XuStand_CylinderLeftFront_log	= new G4LogicalVolume( XuStand_Cylinder,		metalAl,	"XuStand_CylinderLeftFront_log");
	G4LogicalVolume	*XuStand_CylinderLeftBack_log	= new G4LogicalVolume( XuStand_Cylinder,		metalAl,	"XuStand_CylinderLeftBack_log");
	G4LogicalVolume	*XuStand_CylinderRightFront_log	= new G4LogicalVolume( XuStand_Cylinder,		metalAl,	"XuStand_CylinderRightFront_log");
	G4LogicalVolume	*XuStand_CylinderRightBack_log	= new G4LogicalVolume( XuStand_Cylinder,		metalAl,	"XuStand_CylinderRightBack_log");

	G4LogicalVolume	*XuStand_CradlePanelLong_log	= new G4LogicalVolume( XuStand_CradlePanelLong,	metalAl,	"XuStand_CradlePanelLong_log");
	G4LogicalVolume	*XuStand_CradlePanelShort_log	= new G4LogicalVolume( XuStand_CradlePanelShort,metalAl,	"XuStand_CradlePanelShort_log");

	G4LogicalVolume	*XuStand_CradleBarLeftFront_log	= new G4LogicalVolume( XuStand_CradleBarFront,	metalAl,	"XuStand_CradleBarLeftFront_log");
	G4LogicalVolume	*XuStand_CradleBarLeftBack_log	= new G4LogicalVolume( XuStand_CradleBarBack,	metalAl,	"XuStand_CradleBarLeftBack_log");
	G4LogicalVolume	*XuStand_CradleBarRightFront_log= new G4LogicalVolume( XuStand_CradleBarFront,	metalAl,	"XuStand_CradleBarRightFront_log");
	G4LogicalVolume	*XuStand_CradleBarRightBack_log= new G4LogicalVolume( XuStand_CradleBarBack,	metalAl,	"XuStand_CradleBarRightBack_log");

	G4LogicalVolume	*XuStand_Goniometer_log			= new G4LogicalVolume( XuStand_Goniometer,		metalAl,	"XuStand_Goniometer_log");
	// cryostat
	G4LogicalVolume	*OuterCanTube_log				= new G4LogicalVolume( OuterCanTube,			SSteel,		"OuterCanTube_log");
	G4LogicalVolume	*OuterCanBottom_log				= new G4LogicalVolume( OuterCanBottom,			SSteel,		"OuterCanBottom_log");
	G4LogicalVolume	*OuterCanLowerFlange_log		= new G4LogicalVolume( OuterCanLowerFlange,		SSteel,		"OuterCanLowerFlange_log");
	G4LogicalVolume	*OuterCanUpperFlange_log		= new G4LogicalVolume( OuterCanUpperFlange,		SSteel,		"OuterCanUpperFlange_log");

	G4LogicalVolume	*CryostatVacuum_log				= new G4LogicalVolume( CryostatVacuum,			Vacuum,		"CryostatVacuum_log");

	G4LogicalVolume	*AlCanTube_log					= new G4LogicalVolume( AlCanTube,				metalAl, 	"AlCanTube_log");
	G4LogicalVolume	*AlCanBottom_log				= new G4LogicalVolume( AlCanBottom,				metalAl, 	"AlCanBottom_log");

	G4LogicalVolume	*InnerCanTube_log				= new G4LogicalVolume( InnerCanTube,			SSteel, 	"InnerCanTube_log");
	G4LogicalVolume	*InnerCanBottom_log				= new G4LogicalVolume( InnerCanBottom,			SSteel, 	"InnerCanBottom_log");
	G4LogicalVolume	*InnerCanLowerFlange_log		= new G4LogicalVolume( InnerCanLowerFlange,		SSteel, 	"InnerCanLowerFlange_log");
	G4LogicalVolume	*InnerCanUpperFlange_log		= new G4LogicalVolume( InnerCanUpperFlange,		SSteel, 	"InnerCanUpperFlange_log");

	G4LogicalVolume	*SteelHolderUpperSquare_log		= new G4LogicalVolume( SteelHolderUpperSquare1,	SSteel, 	"SteelHolderUpperSquare_log");
	G4LogicalVolume	*SteelHolderLowerSquare_log		= new G4LogicalVolume( SteelHolderLowerSquare1,	SSteel, 	"SteelHolderLowerSquare_log");
	G4LogicalVolume	*SteelHolderTube_log			= new G4LogicalVolume( SteelHolderTube,			SSteel, 	"SteelHolderTube_log");
	// copper finger
	G4LogicalVolume	*CopperFinger_log				= new G4LogicalVolume( CopperFinger,			Copper,		"CopperFinger_log");
	// dewar
	G4LogicalVolume	*DewarOutTube_log				= new G4LogicalVolume( DewarOutTube,			SSteel,		"DewarOutTube_log");
	G4LogicalVolume	*DewarOutBot_log				= new G4LogicalVolume( DewarOutBot,				SSteel,		"DewarOutBot_log");
	G4LogicalVolume	*DewarInTube_log				= new G4LogicalVolume( DewarInTube,				SSteel,		"DewarInTube_log");
	G4LogicalVolume	*DewarInBot_log					= new G4LogicalVolume( DewarInBot,				SSteel,		"DewarInBot_log");
	G4LogicalVolume	*DewarTop_log					= new G4LogicalVolume( DewarTop,				SSteel,		"DewarTop_log");
	// pipes inside the cryostat
	G4LogicalVolume	*Pipe1in_log					= new G4LogicalVolume( Pipe1in,					SSteel,		"Pipe1in_log");
	G4LogicalVolume	*Pipe2in_log					= new G4LogicalVolume( Pipe2in,					SSteel,		"Pipe2in_log");
	G4LogicalVolume	*Pipe3in_log					= new G4LogicalVolume( Pipe3in,					SSteel,		"Pipe3in_log");
	G4LogicalVolume	*Pipe4in_log					= new G4LogicalVolume( Pipe4in,					SSteel,		"Pipe4in_log");
	G4LogicalVolume	*Pipe5in_log					= new G4LogicalVolume( Pipe5in,					SSteel,		"Pipe5in_log");
	G4LogicalVolume	*Pipe6in_log					= new G4LogicalVolume( Pipe6in,					SSteel,		"Pipe6in_log");
	G4LogicalVolume	*Pipe7in_log					= new G4LogicalVolume( Pipe7in,					SSteel,		"Pipe7in_log");
	G4LogicalVolume	*Pipe8in_log					= new G4LogicalVolume( Pipe8in,					SSteel,		"Pipe8in_log");
	G4LogicalVolume	*Pipe9in_log					= new G4LogicalVolume( Pipe9in,					SSteel,		"Pipe9in_log");
	G4LogicalVolume	*Pipe10in_log					= new G4LogicalVolume( Pipe10in,				SSteel,		"Pipe10in_log");
	G4LogicalVolume	*Pipe11in_log					= new G4LogicalVolume( Pipe11in,				SSteel,		"Pipe11in_log");
	// pipes outside the cryostat
	G4LogicalVolume	*Pipe1out_log					= new G4LogicalVolume( Pipe1out,				SSteel,		"Pipe1out_log");
	G4LogicalVolume	*Pipe2out_log					= new G4LogicalVolume( Pipe2out,				SSteel,		"Pipe2out_log");
	G4LogicalVolume	*Pipe3out_log					= new G4LogicalVolume( Pipe3out,				SSteel,		"Pipe3out_log");
	G4LogicalVolume	*Pipe4out_log					= new G4LogicalVolume( Pipe4out,				SSteel,		"Pipe4out_log");
	G4LogicalVolume	*Pipe5out_log					= new G4LogicalVolume( Pipe5out,				SSteel,		"Pipe5out_log");
	G4LogicalVolume	*Pipe6out_log					= new G4LogicalVolume( Pipe6out,				SSteel,		"Pipe6out_log");
	G4LogicalVolume	*Pipe7out_log					= new G4LogicalVolume( Pipe7out,				SSteel,		"Pipe7out_log");
	G4LogicalVolume	*Pipe8out_log					= new G4LogicalVolume( Pipe8out,				SSteel,		"Pipe8out_log");
	G4LogicalVolume	*Pipe9out_log					= new G4LogicalVolume( Pipe9out,				SSteel,		"Pipe9out_log");
	G4LogicalVolume	*Pipe10out_log					= new G4LogicalVolume( Pipe10out,				SSteel,		"Pipe10out_log");
	G4LogicalVolume	*Pipe11out_log					= new G4LogicalVolume( Pipe11out,				SSteel,		"Pipe11out_log");
	G4LogicalVolume	*Pipe1flange_log				= new G4LogicalVolume( Pipe1flange,				SSteel,		"Pipe1flange_log");
	G4LogicalVolume	*Pipe2flange_log				= new G4LogicalVolume( Pipe2flange,				SSteel,		"Pipe2flange_log");
	G4LogicalVolume	*Pipe3flange_log				= new G4LogicalVolume( Pipe3flange,				SSteel,		"Pipe3flange_log");
	G4LogicalVolume	*Pipe4flange_log				= new G4LogicalVolume( Pipe4flange,				SSteel,		"Pipe4flange_log");
	G4LogicalVolume	*Pipe5flange1_log				= new G4LogicalVolume( Pipe5flange,				SSteel,		"Pipe5flange1_log");
	G4LogicalVolume	*Pipe5flange2_log				= new G4LogicalVolume( Pipe5flange,				SSteel,		"Pipe5flange2_log");
	G4LogicalVolume	*Pipe6flange_log				= new G4LogicalVolume( Pipe6flange,				SSteel,		"Pipe6flange_log");
	G4LogicalVolume	*Pipe8flange_log				= new G4LogicalVolume( Pipe8flange,				SSteel,		"Pipe8flange_log");
	G4LogicalVolume	*Pipe9flange_log				= new G4LogicalVolume( Pipe9flange,				SSteel,		"Pipe9flange_log");
	G4LogicalVolume	*Pipe10flange_log				= new G4LogicalVolume( Pipe10flange,			SSteel,		"Pipe10flange_log");
	G4LogicalVolume	*Pipe11flange_log				= new G4LogicalVolume( Pipe11flange,			SSteel,		"Pipe11flange_log");
	// teflon
	G4LogicalVolume	*Support1_log					= new G4LogicalVolume( Support,					Teflon, 	"Support1_log");
	G4LogicalVolume	*Support2_log					= new G4LogicalVolume( Support,					Teflon, 	"Support2_log");
	G4LogicalVolume	*Support3_log					= new G4LogicalVolume( Support,					Teflon, 	"Support3_log");
	G4LogicalVolume	*Support4_log					= new G4LogicalVolume( Support,					Teflon,		"Support4_log");

	G4LogicalVolume *Washer1_log					= new G4LogicalVolume( Washer,					Teflon, 	"Washer1_log");
	G4LogicalVolume	*Washer2_log					= new G4LogicalVolume( Washer,					Teflon, 	"Washer2_log");
	G4LogicalVolume	*Washer3_log					= new G4LogicalVolume( Washer,					Teflon, 	"Washer3_log");
	G4LogicalVolume	*Washer4_log					= new G4LogicalVolume( Washer,					Teflon, 	"Washer4_log");

	G4LogicalVolume *TopClamp_log					= new G4LogicalVolume( TopClamp,				Teflon, 	"TopClamp_log");
	
	G4LogicalVolume	*TopHolderUpperTube_log			= new G4LogicalVolume( TopHolderUpperTube,		Teflon, 	"TopHolderUpperTube_log");
	G4LogicalVolume	*TopHolderLowerTube_log			= new G4LogicalVolume( TopHolderLowerTube,		Teflon, 	"TopHolderLowerTube_log");
	G4LogicalVolume	*TopHolderSquare_log			= new G4LogicalVolume( TopHolderSquare1,		Teflon, 	"TopHolderSquare_log");

	G4LogicalVolume	*ExtractionSpacerUpperHalf_log	= new G4LogicalVolume( ExtractionHalfSpacer,	Teflon, 	"ExtractionSpacerUpperHalf_log");
	G4LogicalVolume *ExtractionSpacerLowerHalf_log	= new G4LogicalVolume( ExtractionHalfSpacer,	Teflon, 	"ExtractionSpacerLowerHalf_log");

	G4LogicalVolume	*DriftSpacerUpperTube_log		= new G4LogicalVolume( DriftSpacerUpperTube,	Teflon, 	"DriftSpacerUpperTube_log");	
	G4LogicalVolume	*DriftSpacerMiddleTube_log		= new G4LogicalVolume( DriftSpacerMiddleTube,	Teflon, 	"DriftSpacerMiddleTube_log");	
	G4LogicalVolume *DriftSpacerLowerTube_log		= new G4LogicalVolume( DriftSpacerLowerTube,	Teflon, 	"DriftSpacerLowerTube_log");	
	G4LogicalVolume	*DriftSpacerUpperSquare_log		= new G4LogicalVolume( DriftSpacerUpperSquare1,	Teflon, 	"DriftSpacerUpperSquare_log");	
	G4LogicalVolume	*DriftSpacerLowerSquare_log		= new G4LogicalVolume( DriftSpacerLowerSquare1,	Teflon, 	"DriftSpacerLowerSquare_log");	

	G4LogicalVolume	*BottomHolderUpperTube_log		= new G4LogicalVolume( BottomHolderUpperTube,	Teflon, 	"BottomHolderUpperTube_log");
	G4LogicalVolume	*BottomHolderMiddleTube_log		= new G4LogicalVolume( BottomHolderMiddleTube,	Teflon, 	"BottomHolderMiddleTube_log");
	G4LogicalVolume	*BottomHolderLowerTube_log		= new G4LogicalVolume( BottomHolderLowerTube,	Teflon, 	"BottomHolderLowerTube_log");

	G4LogicalVolume	*BottomClampUpperPart_log		= new G4LogicalVolume( BottomClampUpperPart,	Teflon, 	"BottomClampUpperPart_log");
	G4LogicalVolume	*BottomClampMiddlePart_log		= new G4LogicalVolume( BottomClampMiddlePart,	Teflon, 	"BottomClampMiddlePart_log");
	G4LogicalVolume	*BottomClampLowerPart_log		= new G4LogicalVolume( BottomClampLowerPart,	Teflon, 	"BottomClampLowerPart_log");

	G4LogicalVolume	*Filler_log						= new G4LogicalVolume( Filler,					Teflon, 	"Filler_log");
	G4LogicalVolume	*topPMTbase_log					= new G4LogicalVolume( PMTbase,					Teflon,		"topPMTbase_log");
	G4LogicalVolume	*bottomPMTbase_log				= new G4LogicalVolume( PMTbase,					Teflon,		"bottomPMTbase_log");

	G4LogicalVolume	*AnodeGrid_log					= new G4LogicalVolume( GridFrame,				SSteel,		"AnodeGrid_log");
	G4LogicalVolume	*CathodeGrid_log				= new G4LogicalVolume( GridFrame,				SSteel,		"CathodeGrid_log");
	G4LogicalVolume	*GateGrid_log					= new G4LogicalVolume( GridFrame,				SSteel,		"GateGrid_log");

	G4LogicalVolume	*topPMTcasing_log				= new G4LogicalVolume( PMTcasing,				SSteel,		"topPMTcasing_log");
	G4LogicalVolume	*topPMTinterior_log				= new G4LogicalVolume( PMTinterior,				Vacuum,		"topPMTinterior_log");
	G4LogicalVolume	*topPMTwindow_log				= new G4LogicalVolume( PMTwindow,				Quartz,		"topPMTwindow_log");
	G4LogicalVolume	*topPMTcathode_log				= new G4LogicalVolume( PMTcathode,				metalAl,	"topPMTcathode_log");

	G4LogicalVolume	*bottomPMTcasing_log			= new G4LogicalVolume( PMTcasing,				SSteel,		"bottomPMTcasing_log");
	G4LogicalVolume *bottomPMTinterior_log			= new G4LogicalVolume( PMTinterior,				Vacuum,		"bottomPMTinterior_log");
	G4LogicalVolume	*bottomPMTwindow_log			= new G4LogicalVolume( PMTwindow,				Quartz,		"bottomPMTwindow_log");
	G4LogicalVolume	*bottomPMTcathode_log			= new G4LogicalVolume( PMTcathode,				metalAl,	"bottomPMTcathode_log");

	// NaI detector, Saint-Gibain mod. 3M3/3
	G4LogicalVolume	*NaI_Crystal_log				= new G4LogicalVolume( NaI_Crystal,				NaI,		"NaI_Crystal_log");
	G4LogicalVolume	*NaI_CrystalHousingTop_log		= new G4LogicalVolume( NaI_CrystalHousingTop,	metalAl,	"NaI_CrystalHousingTop_log");
	G4LogicalVolume	*NaI_CrystalHousingBot_log		= new G4LogicalVolume( NaI_CrystalHousingBot,	metalAl,	"NaI_CrystalHousingBot_log");
	G4LogicalVolume	*NaI_LightShieldTop_log			= new G4LogicalVolume( NaI_LightShieldTop,		MuMetal,	"NaI_LightShieldTop_log");
	G4LogicalVolume	*NaI_LightShieldBot_log			= new G4LogicalVolume( NaI_LightShieldBot,		MuMetal,	"NaI_LightShieldBot_log");
	G4LogicalVolume	*NaI_LightShieldMid_log			= new G4LogicalVolume( NaI_LightShieldMid,		MuMetal,	"NaI_LightShieldMid_log");
	//G4LogicalVolume	*NaI_LightShieldVac_log			= new G4LogicalVolume( NaI_LightShieldVac,		Vacuum,		"NaI_LightShieldVac_log");
	//G4LogicalVolume	*NaI_PMTcasing_log				= new G4LogicalVolume( NaI_PMTcasing,			SSteel,		"NaI_PMTcasing_log");

	// stand for NaI detector
	G4LogicalVolume	*NaIstand_VertBarOut_LeftFront_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"NaIstand_VertBarOut_LeftFront_log");
	G4LogicalVolume	*NaIstand_VertBarOut_LeftBack_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"NaIstand_VertBarOut_LeftBack_log");
	G4LogicalVolume	*NaIstand_VertBarOut_RightFront_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"NaIstand_VertBarOut_RightFront_log");
	G4LogicalVolume	*NaIstand_VertBarOut_RightBack_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"NaIstand_VertBarOut_RightBack_log");

	G4LogicalVolume	*NaIstand_VertBarIn_LeftFront_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"NaIstand_VertBarIn_LeftFront_log");
	G4LogicalVolume	*NaIstand_VertBarIn_LeftBack_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"NaIstand_VertBarIn_LeftBack_log");
	G4LogicalVolume	*NaIstand_VertBarIn_RightFront_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"NaIstand_VertBarIn_RightFront_log");
	G4LogicalVolume	*NaIstand_VertBarIn_RightBack_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"NaIstand_VertBarIn_RightBack_log");

	G4LogicalVolume	*NaIstand_HorzBarOut_Left_log		= new G4LogicalVolume( NaIstand_HorzBarYout,metalAl,	"NaIstand_HorzBarOut_Left_log");
	G4LogicalVolume	*NaIstand_HorzBarOut_Right_log		= new G4LogicalVolume( NaIstand_HorzBarYout,metalAl,	"NaIstand_HorzBarOut_Right_log");
	G4LogicalVolume	*NaIstand_HorzBarOut_Front_log		= new G4LogicalVolume( NaIstand_HorzBarXout,metalAl,	"NaIstand_HorzBarOut_Front_log");
	G4LogicalVolume	*NaIstand_HorzBarOut_Back_log		= new G4LogicalVolume( NaIstand_HorzBarXout,metalAl,	"NaIstand_HorzBarOut_Back_log");

	G4LogicalVolume	*NaIstand_HorzBarIn_Left_log		= new G4LogicalVolume( NaIstand_HorzBarYin,	metalAl,	"NaIstand_HorzBarIn_Left_log");
	G4LogicalVolume	*NaIstand_HorzBarIn_Right_log		= new G4LogicalVolume( NaIstand_HorzBarYin,	metalAl,	"NaIstand_HorzBarIn_Right_log");
	G4LogicalVolume	*NaIstand_HorzBarIn_Front_log		= new G4LogicalVolume( NaIstand_HorzBarXin,	metalAl,	"NaIstand_HorzBarIn_Front_log");
	G4LogicalVolume	*NaIstand_HorzBarIn_Back_log		= new G4LogicalVolume( NaIstand_HorzBarXin,	metalAl,	"NaIstand_HorzBarIn_Back_log");

	G4LogicalVolume	*NaIstand_PlateTop_log				= new G4LogicalVolume( NaIstand_Plate,		metalAl,	"NaIstand_PlateTop_log");
	G4LogicalVolume	*NaIstand_PlateBot_log				= new G4LogicalVolume( NaIstand_Plate,		metalAl,	"NaIstand_PlateBot_log");

	G4LogicalVolume	*NaIstand_CradleBarLeft_log			= new G4LogicalVolume( NaIstand_CradleBar,	metalAl,	"NaIstand_CradleBarLeft_log");
	G4LogicalVolume	*NaIstand_CradleBarRight_log		= new G4LogicalVolume( NaIstand_CradleBar,	metalAl,	"NaIstand_CradleBarRight_log");
	G4LogicalVolume	*NaIstand_CradlePanel_log			= new G4LogicalVolume( NaIstand_CradlePanel,metalAl,	"NaIstand_CradlePanel_log");

	G4LogicalVolume	*NaIcastle_log						= new G4LogicalVolume( NaIcastle,			Lead_m,		"NaIcastle_log");

	// stand for Lead Channel
	G4LogicalVolume	*LeadChannelStand_VertBarOut_LeftFront_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"LeadChannelStand_VertBarOut_LeftFront_log");
	G4LogicalVolume	*LeadChannelStand_VertBarOut_LeftBack_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"LeadChannelStand_VertBarOut_LeftBack_log");
	G4LogicalVolume	*LeadChannelStand_VertBarOut_RightFront_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"LeadChannelStand_VertBarOut_RightFront_log");
	G4LogicalVolume	*LeadChannelStand_VertBarOut_RightBack_log	= new G4LogicalVolume( NaIstand_VertBarOut,	metalAl,	"LeadChannelStand_VertBarOut_RightBack_log");

	G4LogicalVolume	*LeadChannelStand_VertBarIn_LeftFront_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"LeadChannelStand_VertBarIn_LeftFront_log");
	G4LogicalVolume	*LeadChannelStand_VertBarIn_LeftBack_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"LeadChannelStand_VertBarIn_LeftBack_log");
	G4LogicalVolume	*LeadChannelStand_VertBarIn_RightFront_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"LeadChannelStand_VertBarIn_RightFront_log");
	G4LogicalVolume	*LeadChannelStand_VertBarIn_RightBack_log	= new G4LogicalVolume( NaIstand_VertBarIn,	metalAl,	"LeadChannelStand_VertBarIn_RightBack_log");

	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Left_log		= new G4LogicalVolume( NaIstand_HorzBarYout,metalAl,	"LeadChannelStand_HorzBarOut_Left_log");
	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Right_log		= new G4LogicalVolume( NaIstand_HorzBarYout,metalAl,	"LeadChannelStand_HorzBarOut_Right_log");
	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Front_log		= new G4LogicalVolume( NaIstand_HorzBarXout,metalAl,	"LeadChannelStand_HorzBarOut_Front_log");
	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Back_log		= new G4LogicalVolume( NaIstand_HorzBarXout,metalAl,	"LeadChannelStand_HorzBarOut_Back_log");

	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Left_log		= new G4LogicalVolume( NaIstand_HorzBarYin,	metalAl,	"LeadChannelStand_HorzBarIn_Left_log");
	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Right_log		= new G4LogicalVolume( NaIstand_HorzBarYin,	metalAl,	"LeadChannelStand_HorzBarIn_Right_log");
	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Front_log		= new G4LogicalVolume( NaIstand_HorzBarXin,	metalAl,	"LeadChannelStand_HorzBarIn_Front_log");
	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Back_log		= new G4LogicalVolume( NaIstand_HorzBarXin,	metalAl,	"LeadChannelStand_HorzBarIn_Back_log");

	G4LogicalVolume	*LeadChannelStand_PlateTop_log				= new G4LogicalVolume( NaIstand_Plate,		metalAl,	"LeadChannelStand_PlateTop_log");
	G4LogicalVolume	*LeadChannelStand_PlateBot_log				= new G4LogicalVolume( NaIstand_Plate,		metalAl,	"LeadChannelStand_PlateBot_log");
	// collimator for the gamma source
	G4LogicalVolume	*Collimator_log		= new G4LogicalVolume( Collimator,			Lead_m,		"Collimator_log");
	//G4LogicalVolume	*LeadPipe_log	= new G4LogicalVolume( LeadPipe,			Lead_m,		"LeadPipe_log");
	G4LogicalVolume	*LeadChannel_log	= new G4LogicalVolume( LeadChannel,			Lead_m,		"LeadChannel_log");
	G4LogicalVolume	*Aperture_log		= new G4LogicalVolume( Aperture,			Lead_m,		"Aperture_log");
	G4LogicalVolume	*ApertureFrame_log	= new G4LogicalVolume( ApertureFrame,		metalAl,	"ApertureFrame_log");

	//==============================================================================================
	//==============	CONSTRUCTION (Z-dependent from top to bottom)		========================
	//==============================================================================================

	G4double Xuerich_x 		= 0.0 *mm;
	G4double Xuerich_y 		= 0.0 *mm;
	G4double Xuerich_z 		= 273.0 *mm; // top of the upper top flange
	//G4double Xuerich_z 	= 286.5 *mm; // top of the upper top flange
	
	G4double OuterCanUpperFlange_offsetZ	= Xuerich_z - OuterCanUpperFlange_thick/2;
	G4double OuterCanLowerFlange_offsetZ	= OuterCanUpperFlange_offsetZ - OuterCanUpperFlange_thick/2 - OuterCanLowerFlange_thick/2;
	G4double OuterCanTube_offsetZ 			= OuterCanLowerFlange_offsetZ - OuterCanLowerFlange_thick/2 - OuterCanTube_height/2;
	G4double OuterCanBottom_offsetZ			= OuterCanTube_offsetZ - OuterCanTube_height/2 - OuterCanBottom_thick/2;

	G4double PipeIn_offsetZ = OuterCanTube_height/2-InnerPipes_height/2;

	G4double Pipe1_offsetX 			= Xuerich_x - 70 *mm;
	G4double Pipe1_offsetY 			= Xuerich_y + 60 *mm;
	G4double Pipe1out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe1_height/2;
	G4double Pipe1flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe1flange_offset + Pipe1flange_thick/2;

	G4double Pipe2_offsetX 			= Xuerich_x - 30 *mm;
	G4double Pipe2_offsetY 			= Xuerich_y + 55 *mm;
	G4double Pipe2out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe2_height/2;
	G4double Pipe2flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe2flange_offset + Pipe2flange_thick/2;

	G4double Pipe3_offsetX 			= Xuerich_x + 30 *mm;
	G4double Pipe3_offsetY 			= Xuerich_y + 55 *mm;
	G4double Pipe3out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe3_height/2;
	G4double Pipe3flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe3flange_offset + Pipe3flange_thick/2;

	G4double Pipe4_offsetX 			= Xuerich_x + 90 *mm;
	G4double Pipe4_offsetY 			= Xuerich_y + 50 *mm;
	G4double Pipe4out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe4_height/2;
	G4double Pipe4flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe4flange_offset + Pipe4flange_thick/2;
	
	G4double Pipe5_offsetX 			= Xuerich_x - 70 *mm;
	G4double Pipe5_offsetY 			= Xuerich_y + 0.0 *mm;
	G4double Pipe5out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe5_height/2;
	G4double Pipe5flange1_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe5flange1_offset + Pipe5flange_thick/2;
	G4double Pipe5flange2_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe5flange2_offset + Pipe5flange_thick/2;

	G4double Pipe6_offsetX 			= Xuerich_x + 0.0 *mm;
	G4double Pipe6_offsetY 			= Xuerich_y + 0.0 *mm;
	G4double Pipe6out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe6_height/2;
	G4double Pipe6flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe6flange_offset + Pipe6flange_thick/2;

	G4double Pipe7_offsetX 			= Xuerich_x + 70 *mm;
	G4double Pipe7_offsetY 			= Xuerich_y + 0.0 *mm;
	G4double Pipe7out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe7_height/2;

	G4double Pipe8_offsetX 			= Xuerich_x - 30 *mm;
	G4double Pipe8_offsetY 			= Xuerich_y - 55 *mm;
	G4double Pipe8out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe8_height/2;
	G4double Pipe8flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe8flange_offset + Pipe8flange_thick/2;

	G4double Pipe9_offsetX 			= Xuerich_x + 25 *mm;
	G4double Pipe9_offsetY 			= Xuerich_y - 55 *mm;
	G4double Pipe9out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe9_height/2;
	G4double Pipe9flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe9flange_offset + Pipe9flange_thick/2;

	G4double Pipe10_offsetX 		= Xuerich_x + 80 *mm;
	G4double Pipe10_offsetY 		= Xuerich_y - 60 *mm;
	G4double Pipe10out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe10_height/2;
	G4double Pipe10flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe10flange_offset + Pipe10flange_thick/2;

	G4double Pipe11_offsetX 		= Xuerich_x - 75 *mm;
	G4double Pipe11_offsetY 		= Xuerich_y - 60 *mm;
	G4double Pipe11out_offsetZ 		= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe11_height/2;
	G4double Pipe11flange_offsetZ 	= OuterCanUpperFlange_offsetZ + OuterCanUpperFlange_thick/2 + Pipe11flange_offset + Pipe11flange_thick/2;

	G4double PMTpins_length					= 2.0*mm;	
	
	G4double CryostatVacuum_offsetZ			= OuterCanTube_offsetZ;	
	
	G4double InnerCanUpperFlange_offsetZ	= OuterCanTube_height/2 - InnerPipes_height - InnerCanUpperFlange_thick/2;
	G4double InnerCanLowerFlange_offsetZ	= InnerCanUpperFlange_offsetZ - InnerCanLowerFlange_thick/2 - InnerCanLowerFlange_thick/2;
	G4double InnerCanTube_offsetZ			= InnerCanLowerFlange_offsetZ - InnerCanLowerFlange_thick/2 - InnerCanTube_height/2;
	G4double InnerCanBottom_offsetZ			= InnerCanTube_offsetZ - InnerCanTube_height/2 - InnerCanBottom_thick/2;

	// Aluminum can
	G4double AlCanTube_offsetZ 		= InnerCanUpperFlange_offsetZ + InnerCanUpperFlange_thick/2 - AlCanTube_height/2;
	G4double AlCanBottom_offsetZ	= AlCanTube_offsetZ - AlCanTube_height/2 - AlCanBottom_thick/2;
	// copper finger
	G4double CopperFinger_offsetX = Xuerich_x;
	G4double CopperFinger_offsetY = Xuerich_y;
	G4double CopperFinger_offsetZ = OuterCanBottom_offsetZ - OuterCanBottom_thick/2 - CopperFinger_height/2;
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	//--------------	XENON	(SCINTILLATOR CONSTRUCTION - needs to be here, sorry( )	------------////
	G4double GXeVol_height = Washer_thick + SteelHolder_height + Support_height + TopHolderTube_height + GridFrame_thick + ExtractionHalfSpacer_height;
	G4double GXeVol_iR		= 0.0 *mm;
	G4double GXeVol_oR		= InnerCanTube_iR;													
																		
	G4double LXeVol_iR		= 0.0 *mm;															
	G4double LXeVol_oR		= InnerCanTube_iR;													
	G4double LXeVol_height	= InnerCanTube_height + InnerCanLowerFlange_thick - GXeVol_height;	
																							
	G4Tubs *LXeVol		= new G4Tubs("LXeVol", LXeVol_iR, LXeVol_oR, LXeVol_height/2, opendeg, closedeg);																		
	G4Tubs *GXeVol		= new G4Tubs("GXeVol", GXeVol_iR, GXeVol_oR, GXeVol_height/2, opendeg, closedeg);
																				
	G4LogicalVolume	*LXeVol_log		= new G4LogicalVolume( LXeVol,		LXe,	"LXeVol_log");
	G4LogicalVolume	*GXeVol_log		= new G4LogicalVolume( GXeVol,		LXe,	"GXeVol_log");
	////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////
	
	G4double	Washer_offsetZ						= GXeVol_height/2 - Washer_thick/2;
	G4double	Washer1_offsetX						= -35.0 *mm;
	G4double	Washer1_offsetY						=  35.0 *mm;
	G4double	Washer2_offsetX						=  35.0 *mm;
	G4double	Washer2_offsetY						=  35.0 *mm;
	G4double	Washer3_offsetX						= -35.0 *mm;
	G4double	Washer3_offsetY						= -35.0 *mm;
	G4double	Washer4_offsetX						=  35.0 *mm;
	G4double	Washer4_offsetY						= -35.0 *mm;

	G4double	SteelHolderUpperSquare_offsetZ		= Washer_offsetZ - Washer_thick/2 - SteelHolderUpperSquare_thick/2;
	G4double	SteelHolderTube_offsetZ				= Washer_offsetZ - Washer_thick/2 - SteelHolderUpperSquare_thick - SteelHolderTube_height/2;	
	G4double	SteelHolderLowerSquare_offsetZ		= Washer_offsetZ - Washer_thick/2 - SteelHolderUpperSquare_thick - SteelHolderTube_height - SteelHolderLowerSquare_thick/2;		
	
	G4double	Support_offsetZ						= SteelHolderLowerSquare_offsetZ - SteelHolderLowerSquare_thick/2 - Support_height/2;
	G4double	Support1_offsetX					= -35.5 *mm;
	G4double	Support1_offsetY					=  35.5 *mm;
	G4double	Support2_offsetX					=  35.5 *mm;
	G4double	Support2_offsetY					=  35.5 *mm;
	G4double	Support3_offsetX					= -35.5 *mm;
	G4double	Support3_offsetY					= -35.5 *mm;
	G4double	Support4_offsetX					=  35.5 *mm;
	G4double	Support4_offsetY					= -35.5 *mm;

	G4double	TopHolderUpperTube_offsetZ			= Support_offsetZ - Support_height/2 - TopHolderUpperTube_height/2;
	G4double	TopHolderLowerTube_offsetZ			= TopHolderUpperTube_offsetZ - TopHolderUpperTube_height/2 - TopHolderLowerTube_height/2;
	G4double	TopHolderSquare_offsetZ				= TopHolderLowerTube_offsetZ - TopHolderLowerTube_height/2 - TopHolderSquare_thick/2;
	
	G4double	GapForWires							= 0.5 *mm;
	G4double	AnodeGrid_offsetZ					= TopHolderLowerTube_offsetZ - TopHolderLowerTube_height/2 - GridFrame_thick/2;
	G4double	ExtractionSpacerUpperHalf_offsetZ	= AnodeGrid_offsetZ - GridFrame_thick/2 - ExtractionHalfSpacer_height/2;
		
	G4double	ExtractionSpacerLowerHalf_offsetZ	= LXeVol_height/2 - ExtractionHalfSpacer_height/2;
	G4double	GateGrid_offsetZ					= ExtractionSpacerLowerHalf_offsetZ - ExtractionHalfSpacer_height/2 - GridFrame_thick/2;
	
	G4double	DriftSpacerUpperSquare_offsetZ		= GateGrid_offsetZ - GridFrame_thick/2 + DriftSpacerUpperSquare_thick/2;
	G4double	DriftSpacerUpperTube_offsetZ		= DriftSpacerUpperSquare_offsetZ - DriftSpacerUpperSquare_thick/2 - DriftSpacerUpperTube_height/2;
	G4double	DriftSpacerMiddleTube_offsetZ		= DriftSpacerUpperTube_offsetZ - DriftSpacerUpperTube_height/2 - DriftSpacerMiddleTube_height/2;
	G4double	DriftSpacerLowerTube_offsetZ		= DriftSpacerMiddleTube_offsetZ - DriftSpacerMiddleTube_height/2 - DriftSpacerLowerTube_height/2;
	G4double	DriftSpacerLowerSquare_offsetZ		= DriftSpacerLowerTube_offsetZ - DriftSpacerLowerTube_height/2 - DriftSpacerLowerSquare_thick/2;
	
	G4double	CathodeGrid_offsetZ					= DriftSpacerLowerTube_offsetZ - DriftSpacerLowerTube_height/2 - GridFrame_thick/2;

	G4double	BottomHolderUpperTube_offsetZ		= CathodeGrid_offsetZ - GridFrame_thick/2 - BottomHolderUpperTube_height/2;
	G4double	BottomHolderMiddleTube_offsetZ		= BottomHolderUpperTube_offsetZ - BottomHolderUpperTube_height/2 - BottomHolderMiddleTube_height/2;
	G4double	BottomHolderLowerTube_offsetZ		= BottomHolderMiddleTube_offsetZ - BottomHolderMiddleTube_height/2 - BottomHolderLowerTube_height/2;
	
	G4double	Filler_offsetZ						= BottomHolderLowerTube_offsetZ - BottomHolderLowerTube_height/2 + Filler_height/2;
	
	G4double	topPMTcasing_offsetZ				= TopHolderUpperTube_offsetZ - TopHolderUpperTube_height/2 + PMTcasing_height/2;
	G4double	topPMTinterior_offsetZ				= -PMTcasing_height/2 + PMTinterior_height/2;
	G4double	topPMTwindow_offsetZ				= -PMTinterior_height/2 + PMTwindow_thick/2;
	G4double	topPMTcathode_offsetZ				= PMTwindow_thick/2 - PMTcathode_thick/2;
	
	G4double	topPMTbase_offsetZ					= topPMTcasing_offsetZ + PMTcasing_height/2 + PMTpins_length + PMTbase_thick/2;	
	G4double	TopClamp_offsetZ					= topPMTcasing_offsetZ + PMTcasing_height/2 + TopClamp_thick/2;

	G4double	bottomPMTcasing_offsetZ				= BottomHolderMiddleTube_offsetZ + BottomHolderMiddleTube_height/2 - PMTcasing_height/2;
	G4double	bottomPMTinterior_offsetZ			= PMTcasing_height/2 - PMTinterior_height/2;
	G4double	bottomPMTwindow_offsetZ				= PMTinterior_height/2 - PMTwindow_thick/2;
	G4double	bottomPMTcathode_offsetZ			= -PMTwindow_thick/2 + PMTcathode_thick/2;

	G4double	BottomClampUpperPart_offsetZ		= bottomPMTcasing_offsetZ - PMTcasing_height/2 - BottomClampUpperPart_height/2;
	G4double	BottomClampMiddlePart_offsetZ		= BottomClampUpperPart_offsetZ - BottomClampUpperPart_height/2 - BottomClampMiddlePart_height/2;
	G4double	BottomClampLowerPart_offsetZ		= BottomClampMiddlePart_offsetZ - BottomClampMiddlePart_height/2 - BottomClampLowerPart_height/2;
	
	G4double	bottomPMTbase_offsetZ				= bottomPMTcasing_offsetZ - PMTcasing_height/2 - PMTpins_length - PMTbase_thick/2;
	
	/////////////////////////////////////////////////////////////////////////////////////////////
	// TARGET VOLUME
	G4double LXeTarget_iR		= 0.0 *mm;
	G4double LXeTarget_oR		= DriftSpacerMiddleTube_iR;
	G4double LXeTarget_height	= DriftSpacerLowerTube_height + DriftSpacerMiddleTube_height + DriftSpacerUpperTube_height + ExtractionHalfSpacer_height + GridFrame_thick*2 + BottomHolderUpperTube_height;	

	G4double LXeGate_iR			= 0.0 *mm;															
	G4double LXeGate_oR			= ExtractionSpacer_iR;	
	G4double LXeGate_height		= ExtractionHalfSpacer_height + GridFrame_thick + TopHolderLowerTube_height;	

	//G4double LXeCathode_iR		= 0.0 *mm;															
	//G4double LXeCathode_oR		= BottomHolderUpperTube_iR;;													
	//G4double LXeCathode_height	= BottomHolderUpperTube_height;	

	G4double LXeVol_offsetZ	= InnerCanTube_offsetZ - InnerCanTube_height/2 + LXeVol_height/2;	
	G4double GXeVol_offsetZ	= InnerCanTube_offsetZ + InnerCanLowerFlange_thick + InnerCanTube_height/2 - GXeVol_height/2;

	//G4double LXeTarget_offsetZ	= DriftSpacerUpperSquare_offsetZ + DriftSpacerUpperSquare_thick/2 - LXeTarget_height/2;
	G4double LXeTarget_offsetZ		= LXeVol_height/2 - LXeTarget_height/2;
	//G4double LXeGate_offsetZ		= ExtractionSpacerLowerHalf_offsetZ;
	G4double LXeGate_offsetZ		= -GXeVol_height/2 + LXeGate_height/2;
	//G4double LXeCathode_offsetZ	= BottomHolderUpperTube_offsetZ;
	
	G4Tubs *LXeTarget	= new G4Tubs("LXeTarget", 	LXeTarget_iR, LXeTarget_oR, LXeTarget_height/2, opendeg, closedeg);																		
	G4Tubs *LXeGate		= new G4Tubs("LXeGate", 	LXeGate_iR, LXeGate_oR, LXeGate_height/2, opendeg, closedeg);																		
	//G4Tubs *LXeCathode	= new G4Tubs("LXeCathode", 	LXeCathode_iR, LXeCathode_oR, LXeCathode_height/2, opendeg, closedeg);																		

	G4LogicalVolume	*LXeTarget_log	= new G4LogicalVolume( LXeTarget,	LXe,	"LXeTarget_log");
	G4LogicalVolume	*LXeGate_log	= new G4LogicalVolume( LXeGate,		LXe,	"LXeGate_log");
	//G4LogicalVolume	*LXeCathode_log	= new G4LogicalVolume( LXeCathode,	LXe,	"LXeCathode_log");
	///////////////////////////////////////////////////////////////////////////////////////////////


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// DEFINE ROTATION ANGLES
	G4double NaIcastle_offsetDepth 	= 8. *cm;
	//G4double NaI_distance 		= 108. *cm + OuterCanTube_oR - ApertureFrame_thick;
	G4double NaI_distance 			= XuStand_width/2 + LeadChannel_length + NaIcastle_offsetDepth + NaI_Crystal_height/2;

	G4double ScatterAngle_deg 		= -0.0;		// *deg;
	G4double ScatterAngle_degA 		= 0.0 *deg;	// *deg;
	G4double ScatterAngle_rad 		= ScatterAngle_deg * TMath::Pi() / 180;// *rad;
	
	G4RotationMatrix RotationNaI;
	RotationNaI.rotateY(-ScatterAngle_degA);
	RotationNaI.rotateX(90.0 *deg);
	RotationNaI.rotateZ(0.0 *deg);

	G4RotationMatrix RotationNaIcastle;
	RotationNaIcastle.rotateZ(-ScatterAngle_degA);
	RotationNaIcastle.rotateX(0.0 *deg);
	RotationNaIcastle.rotateY(0.0 *deg);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// XUERICH DEFAULT STAND, ROTATED
	G4double XuStand_HorzBar_distance = XuStand_width/2 - XuStand_profile_out/2;
	
	G4double XuStand_HorzBar_TopLeft_offsetX 	= Xuerich_x - XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopLeft_offsetY 	= Xuerich_y - XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopLeft_offsetZ 	= Xuerich_z - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_TopLeft_V (XuStand_HorzBar_TopLeft_offsetX, XuStand_HorzBar_TopLeft_offsetY, XuStand_HorzBar_TopLeft_offsetZ);
	G4Transform3D XuStand_HorzBar_TopLeft_T(RotationNaIcastle, XuStand_HorzBar_TopLeft_V);	

	G4double XuStand_HorzBar_TopRight_offsetX 	= Xuerich_x + XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopRight_offsetY 	= Xuerich_y + XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopRight_offsetZ 	= Xuerich_z - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_TopRight_V (XuStand_HorzBar_TopRight_offsetX, XuStand_HorzBar_TopRight_offsetY, XuStand_HorzBar_TopRight_offsetZ);
	G4Transform3D XuStand_HorzBar_TopRight_T(RotationNaIcastle, XuStand_HorzBar_TopRight_V);	

	G4double XuStand_HorzBar_TopFront_offsetX 	= Xuerich_x + XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopFront_offsetY 	= Xuerich_y - XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopFront_offsetZ 	= Xuerich_z - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_TopFront_V (XuStand_HorzBar_TopFront_offsetX, XuStand_HorzBar_TopFront_offsetY, XuStand_HorzBar_TopFront_offsetZ);
	G4Transform3D XuStand_HorzBar_TopFront_T(RotationNaIcastle, XuStand_HorzBar_TopFront_V);	

	G4double XuStand_HorzBar_TopBack_offsetX 	= Xuerich_x - XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopBack_offsetY 	= Xuerich_y + XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_TopBack_offsetZ 	= Xuerich_z - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_TopBack_V (XuStand_HorzBar_TopBack_offsetX, XuStand_HorzBar_TopBack_offsetY, XuStand_HorzBar_TopBack_offsetZ);
	G4Transform3D XuStand_HorzBar_TopBack_T(RotationNaIcastle, XuStand_HorzBar_TopBack_V);

	G4double XuStand_Wing_distance = XuStand_depth/2 - XuStand_Wing_length/2;

	G4double XuStand_WingFront_offsetX 			= Xuerich_x + XuStand_Wing_distance*sin(ScatterAngle_rad);
	G4double XuStand_WingFront_offsetY 			= Xuerich_y - XuStand_Wing_distance*cos(ScatterAngle_rad);
	G4double XuStand_WingFront_offsetZ 			= Xuerich_z + XuStand_Wing_thick/2;
	G4ThreeVector XuStand_WingFront_V (XuStand_WingFront_offsetX, XuStand_WingFront_offsetY, XuStand_WingFront_offsetZ);
	G4Transform3D XuStand_WingFront_T(RotationNaIcastle, XuStand_WingFront_V);

	G4double XuStand_WingBack_offsetX 			= Xuerich_x - XuStand_Wing_distance*sin(ScatterAngle_rad);
	G4double XuStand_WingBack_offsetY 			= Xuerich_y + XuStand_Wing_distance*cos(ScatterAngle_rad);
	G4double XuStand_WingBack_offsetZ 			= Xuerich_z + XuStand_Wing_thick/2;
	G4ThreeVector XuStand_WingBack_V (XuStand_WingBack_offsetX, XuStand_WingBack_offsetY, XuStand_WingBack_offsetZ);
	G4Transform3D XuStand_WingBack_T(RotationNaIcastle, XuStand_WingBack_V);

	G4double XuStand_HorzBar_MidLeft_offsetX 	= Xuerich_x - XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidLeft_offsetY 	= Xuerich_y - XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidLeft_offsetZ 	= Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_MidLeft_V (XuStand_HorzBar_MidLeft_offsetX, XuStand_HorzBar_MidLeft_offsetY, XuStand_HorzBar_MidLeft_offsetZ);
	G4Transform3D XuStand_HorzBar_MidLeft_T(RotationNaIcastle, XuStand_HorzBar_MidLeft_V);

	G4double XuStand_HorzBar_MidRight_offsetX 	= Xuerich_x + XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidRight_offsetY 	= Xuerich_y + XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidRight_offsetZ 	= Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_MidRight_V (XuStand_HorzBar_MidRight_offsetX, XuStand_HorzBar_MidRight_offsetY, XuStand_HorzBar_MidRight_offsetZ);
	G4Transform3D XuStand_HorzBar_MidRight_T(RotationNaIcastle, XuStand_HorzBar_MidRight_V);

	G4double XuStand_HorzBar_MidFront_offsetX 	= Xuerich_x + XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidFront_offsetY 	= Xuerich_y - XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidFront_offsetZ 	= Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_MidFront_V (XuStand_HorzBar_MidFront_offsetX, XuStand_HorzBar_MidFront_offsetY, XuStand_HorzBar_MidFront_offsetZ);
	G4Transform3D XuStand_HorzBar_MidFront_T(RotationNaIcastle, XuStand_HorzBar_MidFront_V);

	G4double XuStand_HorzBar_MidBack_offsetX 	= Xuerich_x - XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidBack_offsetY 	= Xuerich_y + XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_MidBack_offsetZ 	= Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_MidBack_V (XuStand_HorzBar_MidBack_offsetX, XuStand_HorzBar_MidBack_offsetY, XuStand_HorzBar_MidBack_offsetZ);
	G4Transform3D XuStand_HorzBar_MidBack_T(RotationNaIcastle, XuStand_HorzBar_MidBack_V);

	G4double XuStand_VertBar_LeftFront_offsetX 	= Xuerich_x - (XuStand_HorzBar_distance*cos(ScatterAngle_rad)-XuStand_HorzBar_distance*sin(ScatterAngle_rad));
	G4double XuStand_VertBar_LeftFront_offsetY 	= Xuerich_y - (XuStand_HorzBar_distance*sin(ScatterAngle_rad)+XuStand_HorzBar_distance*cos(ScatterAngle_rad));
	G4double XuStand_VertBar_LeftFront_offsetZ 	= Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out - XuStand_VertBar_height/2;
	G4ThreeVector XuStand_VertBar_LeftFront_V (XuStand_VertBar_LeftFront_offsetX, XuStand_VertBar_LeftFront_offsetY, XuStand_VertBar_LeftFront_offsetZ);
	G4Transform3D XuStand_VertBar_LeftFront_T(RotationNaIcastle, XuStand_VertBar_LeftFront_V);

	G4double XuStand_VertBar_LeftBack_offsetX 	= Xuerich_x - (XuStand_HorzBar_distance*cos(ScatterAngle_rad)+XuStand_HorzBar_distance*sin(ScatterAngle_rad));
	G4double XuStand_VertBar_LeftBack_offsetY 	= Xuerich_y - (XuStand_HorzBar_distance*sin(ScatterAngle_rad)-XuStand_HorzBar_distance*cos(ScatterAngle_rad));
	G4double XuStand_VertBar_LeftBack_offsetZ 	= Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out - XuStand_VertBar_height/2;
	G4ThreeVector XuStand_VertBar_LeftBack_V (XuStand_VertBar_LeftBack_offsetX, XuStand_VertBar_LeftBack_offsetY, XuStand_VertBar_LeftBack_offsetZ);
	G4Transform3D XuStand_VertBar_LeftBack_T(RotationNaIcastle, XuStand_VertBar_LeftBack_V);

	G4double XuStand_VertBar_RightFront_offsetX = Xuerich_x + (XuStand_HorzBar_distance*cos(ScatterAngle_rad)+XuStand_HorzBar_distance*sin(ScatterAngle_rad));
	G4double XuStand_VertBar_RightFront_offsetY = Xuerich_y + (XuStand_HorzBar_distance*sin(ScatterAngle_rad)-XuStand_HorzBar_distance*cos(ScatterAngle_rad));
	G4double XuStand_VertBar_RightFront_offsetZ = Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out - XuStand_VertBar_height/2;
	G4ThreeVector XuStand_VertBar_RightFront_V (XuStand_VertBar_RightFront_offsetX, XuStand_VertBar_RightFront_offsetY, XuStand_VertBar_RightFront_offsetZ);
	G4Transform3D XuStand_VertBar_RightFront_T(RotationNaIcastle, XuStand_VertBar_RightFront_V);

	G4double XuStand_VertBar_RightBack_offsetX 	= Xuerich_x + (XuStand_HorzBar_distance*cos(ScatterAngle_rad)-XuStand_HorzBar_distance*sin(ScatterAngle_rad));
	G4double XuStand_VertBar_RightBack_offsetY 	= Xuerich_y + (XuStand_HorzBar_distance*sin(ScatterAngle_rad)+XuStand_HorzBar_distance*cos(ScatterAngle_rad));
	G4double XuStand_VertBar_RightBack_offsetZ 	= Xuerich_z - XuStand_profile_out - XuStand_GapTopMid - XuStand_profile_out - XuStand_VertBar_height/2;
	G4ThreeVector XuStand_VertBar_RightBack_V (XuStand_VertBar_RightBack_offsetX, XuStand_VertBar_RightBack_offsetY, XuStand_VertBar_RightBack_offsetZ);
	G4Transform3D XuStand_VertBar_RightBack_T(RotationNaIcastle, XuStand_VertBar_RightBack_V);

	G4double XuStand_BottomPlate_offsetX 		= Xuerich_x;
	G4double XuStand_BottomPlate_offsetY 		= Xuerich_y;
	G4double XuStand_BottomPlate_offsetZ 		= XuStand_VertBar_LeftFront_offsetZ  - XuStand_VertBar_height/2 - XuStand_BottomPlate_thick/2;
	G4ThreeVector XuStand_BottomPlate_V (XuStand_BottomPlate_offsetX, XuStand_BottomPlate_offsetY, XuStand_BottomPlate_offsetZ);
	G4Transform3D XuStand_BottomPlate_T(RotationNaIcastle, XuStand_BottomPlate_V);

	G4double XuStand_HorzBar_BotLeft_offsetX 	= Xuerich_x - XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotLeft_offsetY 	= Xuerich_y - XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotLeft_offsetZ 	= XuStand_BottomPlate_offsetZ + XuStand_BottomPlate_thick/2 + XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_BotLeft_V (XuStand_HorzBar_BotLeft_offsetX, XuStand_HorzBar_BotLeft_offsetY, XuStand_HorzBar_BotLeft_offsetZ);
	G4Transform3D XuStand_HorzBar_BotLeft_T(RotationNaIcastle, XuStand_HorzBar_BotLeft_V);

	G4double XuStand_HorzBar_BotRight_offsetX 	= Xuerich_x + XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotRight_offsetY 	= Xuerich_y + XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotRight_offsetZ 	= XuStand_BottomPlate_offsetZ + XuStand_BottomPlate_thick/2 + XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_BotRight_V (XuStand_HorzBar_BotRight_offsetX, XuStand_HorzBar_BotRight_offsetY, XuStand_HorzBar_BotRight_offsetZ);
	G4Transform3D XuStand_HorzBar_BotRight_T(RotationNaIcastle, XuStand_HorzBar_BotRight_V);

	G4double XuStand_HorzBar_BotFront_offsetX 	= Xuerich_x + XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotFront_offsetY 	= Xuerich_y - XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotFront_offsetZ 	= XuStand_BottomPlate_offsetZ + XuStand_BottomPlate_thick/2 + XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_BotFront_V (XuStand_HorzBar_BotFront_offsetX, XuStand_HorzBar_BotFront_offsetY, XuStand_HorzBar_BotFront_offsetZ);
	G4Transform3D XuStand_HorzBar_BotFront_T(RotationNaIcastle, XuStand_HorzBar_BotFront_V);

	G4double XuStand_HorzBar_BotBack_offsetX 	= Xuerich_x - XuStand_HorzBar_distance*sin(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotBack_offsetY 	= Xuerich_y + XuStand_HorzBar_distance*cos(ScatterAngle_rad);
	G4double XuStand_HorzBar_BotBack_offsetZ 	= XuStand_BottomPlate_offsetZ + XuStand_BottomPlate_thick/2 + XuStand_profile_out/2;
	G4ThreeVector XuStand_HorzBar_BotBack_V (XuStand_HorzBar_BotBack_offsetX, XuStand_HorzBar_BotBack_offsetY, XuStand_HorzBar_BotBack_offsetZ);
	G4Transform3D XuStand_HorzBar_BotBack_T(RotationNaIcastle, XuStand_HorzBar_BotBack_V);

	G4double XuStand_Cylinder_distanceX 		= XuStand_HorzBar_distance;
	G4double XuStand_Cylinder_distanceY 		= XuStand_CradlePanelLong_length/2 - XuStand_Cylinder_oD/2;

	G4double XuStand_CylinderLeftFront_offsetX 	= Xuerich_x - (XuStand_Cylinder_distanceX*cos(ScatterAngle_rad)-XuStand_Cylinder_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CylinderLeftFront_offsetY 	= Xuerich_y - (XuStand_Cylinder_distanceX*sin(ScatterAngle_rad)+XuStand_Cylinder_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CylinderLeftFront_offsetZ 	= XuStand_HorzBar_TopLeft_offsetZ + XuStand_profile_out/2 + XuStand_Cylinder_height/2;
	G4ThreeVector XuStand_CylinderLeftFront_V (XuStand_CylinderLeftFront_offsetX, XuStand_CylinderLeftFront_offsetY, XuStand_CylinderLeftFront_offsetZ);
	G4Transform3D XuStand_CylinderLeftFront_T(RotationNaIcastle, XuStand_CylinderLeftFront_V);

	G4double XuStand_CylinderLeftBack_offsetX 	= Xuerich_x - (XuStand_Cylinder_distanceX*cos(ScatterAngle_rad)+XuStand_Cylinder_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CylinderLeftBack_offsetY 	= Xuerich_y - (XuStand_Cylinder_distanceX*sin(ScatterAngle_rad)-XuStand_Cylinder_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CylinderLeftBack_offsetZ 	= XuStand_HorzBar_TopLeft_offsetZ + XuStand_profile_out/2 + XuStand_Cylinder_height/2;
	G4ThreeVector XuStand_CylinderLeftBack_V (XuStand_CylinderLeftBack_offsetX, XuStand_CylinderLeftBack_offsetY, XuStand_CylinderLeftBack_offsetZ);
	G4Transform3D XuStand_CylinderLeftBack_T(RotationNaIcastle, XuStand_CylinderLeftBack_V);

	G4double XuStand_CylinderRightFront_offsetX = Xuerich_x + (XuStand_Cylinder_distanceX*cos(ScatterAngle_rad)+XuStand_Cylinder_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CylinderRightFront_offsetY = Xuerich_y + (XuStand_Cylinder_distanceX*sin(ScatterAngle_rad)-XuStand_Cylinder_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CylinderRightFront_offsetZ = XuStand_HorzBar_TopRight_offsetZ + XuStand_profile_out/2 + XuStand_Cylinder_height/2;
	G4ThreeVector XuStand_CylinderRightFront_V (XuStand_CylinderRightFront_offsetX, XuStand_CylinderRightFront_offsetY, XuStand_CylinderRightFront_offsetZ);
	G4Transform3D XuStand_CylinderRightFront_T(RotationNaIcastle, XuStand_CylinderRightFront_V);

	G4double XuStand_CylinderRightBack_offsetX 	= Xuerich_x + (XuStand_Cylinder_distanceX*cos(ScatterAngle_rad)-XuStand_Cylinder_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CylinderRightBack_offsetY 	= Xuerich_y + (XuStand_Cylinder_distanceX*sin(ScatterAngle_rad)+XuStand_Cylinder_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CylinderRightBack_offsetZ 	= XuStand_HorzBar_TopRight_offsetZ + XuStand_profile_out/2 + XuStand_Cylinder_height/2;
	G4ThreeVector XuStand_CylinderRightBack_V (XuStand_CylinderRightBack_offsetX, XuStand_CylinderRightBack_offsetY, XuStand_CylinderRightBack_offsetZ);
	G4Transform3D XuStand_CylinderRightBack_T(RotationNaIcastle, XuStand_CylinderRightBack_V);

	G4double XuStand_CradleBar_distanceX 		= XuStand_HorzBar_distance;
	G4double XuStand_CradleBar_distanceY 		= XuStand_CradlePanelLong_length/2 - XuStand_CradleBarBot_length + XuStand_CradleBarTop_oD/2;

	G4double XuStand_CradleBarLeftFront_offsetX = Xuerich_x - (XuStand_CradleBar_distanceX*cos(ScatterAngle_rad)-XuStand_CradleBar_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CradleBarLeftFront_offsetY = Xuerich_y - (XuStand_CradleBar_distanceX*sin(ScatterAngle_rad)+XuStand_CradleBar_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CradleBarLeftFront_offsetZ = XuStand_CylinderLeftFront_offsetZ + XuStand_Cylinder_height/2 + XuStand_CradleBarBot_thick + XuStand_CradleBarTop_height/2;
	G4ThreeVector XuStand_CradleBarLeftFront_V (XuStand_CradleBarLeftFront_offsetX, XuStand_CradleBarLeftFront_offsetY, XuStand_CradleBarLeftFront_offsetZ);
	G4Transform3D XuStand_CradleBarLeftFront_T(RotationNaIcastle, XuStand_CradleBarLeftFront_V);

	G4double XuStand_CradleBarLeftBack_offsetX 	= Xuerich_x - (XuStand_CradleBar_distanceX*cos(ScatterAngle_rad)+XuStand_CradleBar_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CradleBarLeftBack_offsetY 	= Xuerich_y - (XuStand_CradleBar_distanceX*sin(ScatterAngle_rad)-XuStand_CradleBar_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CradleBarLeftBack_offsetZ 	= XuStand_CylinderLeftBack_offsetZ + XuStand_Cylinder_height/2 + XuStand_CradleBarBot_thick + XuStand_CradleBarTop_height/2;
	G4ThreeVector XuStand_CradleBarLeftBack_V (XuStand_CradleBarLeftBack_offsetX, XuStand_CradleBarLeftBack_offsetY, XuStand_CradleBarLeftBack_offsetZ);
	G4Transform3D XuStand_CradleBarLeftBack_T(RotationNaIcastle, XuStand_CradleBarLeftBack_V);

	G4double XuStand_CradleBarRightFront_offsetX = Xuerich_x + (XuStand_CradleBar_distanceX*cos(ScatterAngle_rad)+XuStand_CradleBar_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CradleBarRightFront_offsetY = Xuerich_y + (XuStand_CradleBar_distanceX*sin(ScatterAngle_rad)-XuStand_CradleBar_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CradleBarRightFront_offsetZ = XuStand_CylinderRightFront_offsetZ + XuStand_Cylinder_height/2 + XuStand_CradleBarBot_thick + XuStand_CradleBarTop_height/2;
	G4ThreeVector XuStand_CradleBarRightFront_V (XuStand_CradleBarRightFront_offsetX, XuStand_CradleBarRightFront_offsetY, XuStand_CradleBarRightFront_offsetZ);
	G4Transform3D XuStand_CradleBarRightFront_T(RotationNaIcastle, XuStand_CradleBarRightFront_V);

	G4double XuStand_CradleBarRightBack_offsetX = Xuerich_x + (XuStand_CradleBar_distanceX*cos(ScatterAngle_rad)-XuStand_CradleBar_distanceY*sin(ScatterAngle_rad));
	G4double XuStand_CradleBarRightBack_offsetY = Xuerich_y + (XuStand_CradleBar_distanceX*sin(ScatterAngle_rad)+XuStand_CradleBar_distanceY*cos(ScatterAngle_rad));
	G4double XuStand_CradleBarRightBack_offsetZ = XuStand_CylinderRightBack_offsetZ + XuStand_Cylinder_height/2 + XuStand_CradleBarBot_thick + XuStand_CradleBarTop_height/2;
	G4ThreeVector XuStand_CradleBarRightBack_V (XuStand_CradleBarRightBack_offsetX, XuStand_CradleBarRightBack_offsetY, XuStand_CradleBarRightBack_offsetZ);
	G4Transform3D XuStand_CradleBarRightBack_T(RotationNaIcastle, XuStand_CradleBarRightBack_V);

	G4double XuStand_CradlePanel_distance = XuStand_HorzBar_distance;

	G4double XuStand_CradlePanelLong_offsetX 	= Xuerich_x - XuStand_CradlePanel_distance*cos(ScatterAngle_rad);
	G4double XuStand_CradlePanelLong_offsetY 	= Xuerich_y - XuStand_CradlePanel_distance*sin(ScatterAngle_rad);
	G4double XuStand_CradlePanelLong_offsetZ 	= XuStand_CylinderLeftFront_offsetZ + XuStand_Cylinder_height/2 + XuStand_CradleBar_height + XuStand_CradlePanel_thick/2;
	G4ThreeVector XuStand_CradlePanelLong_V (XuStand_CradlePanelLong_offsetX, XuStand_CradlePanelLong_offsetY, XuStand_CradlePanelLong_offsetZ);
	G4Transform3D XuStand_CradlePanelLong_T(RotationNaIcastle, XuStand_CradlePanelLong_V);

	G4double XuStand_CradlePanelShort_offsetX 	= Xuerich_x + XuStand_CradlePanel_distance*cos(ScatterAngle_rad);
	G4double XuStand_CradlePanelShort_offsetY 	= Xuerich_y + XuStand_CradlePanel_distance*sin(ScatterAngle_rad);
	G4double XuStand_CradlePanelShort_offsetZ 	= XuStand_CylinderLeftFront_offsetZ + XuStand_Cylinder_height/2 + XuStand_CradleBar_height + XuStand_CradlePanel_thick/2;
	G4ThreeVector XuStand_CradlePanelShort_V (XuStand_CradlePanelShort_offsetX, XuStand_CradlePanelShort_offsetY, XuStand_CradlePanelShort_offsetZ);
	G4Transform3D XuStand_CradlePanelShort_T(RotationNaIcastle, XuStand_CradlePanelShort_V);

	G4double XuStand_Goniometer_offsetX 		= Xuerich_x;
	G4double XuStand_Goniometer_offsetY 		= Xuerich_y;
	G4double XuStand_Goniometer_offsetZ 		= XuStand_CradlePanelLong_offsetZ + XuStand_CradlePanel_thick + XuStand_Goniometer_thick/2;

	// dewar
	G4double DewarOutBot_offsetX 	= Xuerich_x;
	G4double DewarOutBot_offsetY 	= Xuerich_y;
	G4double DewarOutBot_offsetZ 	= XuStand_BottomPlate_offsetZ + XuStand_BottomPlate_thick/2 + Dewar_thick/2;

	G4double DewarOutTube_offsetX 	= Xuerich_x;
	G4double DewarOutTube_offsetY 	= Xuerich_y;
	G4double DewarOutTube_offsetZ 	= DewarOutBot_offsetZ + Dewar_thick/2 + DewarOutTube_height/2;

	G4double DewarInTube_offsetX 	= Xuerich_x;
	G4double DewarInTube_offsetY 	= Xuerich_y;
	G4double DewarInTube_offsetZ 	= DewarOutTube_offsetZ + DewarOutTube_height/2 - DewarInTube_height/2;

	G4double DewarInBot_offsetX 	= Xuerich_x;
	G4double DewarInBot_offsetY 	= Xuerich_y;
	G4double DewarInBot_offsetZ 	= DewarInTube_offsetZ - DewarOutTube_height/2 - Dewar_thick/2;

	G4double DewarTop_offsetX 		= Xuerich_x;
	G4double DewarTop_offsetY 		= Xuerich_y;
	G4double DewarTop_offsetZ 		= DewarOutTube_offsetZ + DewarOutTube_height/2 + Dewar_thick/2;

	// NaI detector, Saint-Gibain mod. 3M3/3
	G4double NaI_x	= Xuerich_x*cos(ScatterAngle_rad)-(Xuerich_y-NaI_distance)*sin(ScatterAngle_rad);
	G4double NaI_y	= Xuerich_x*sin(ScatterAngle_rad)+(Xuerich_y-NaI_distance)*cos(ScatterAngle_rad);
	G4double NaI_z	= 0.0 *mm;
	
	// horizontal arrangement, facing Xuerich, related to the center of the NaI crystal
	G4double NaI_Crystal_offsetX = NaI_x;
	G4double NaI_Crystal_offsetY = NaI_y;
	G4double NaI_Crystal_offsetZ = NaI_z;
	G4ThreeVector NaI_Crystal_V (NaI_Crystal_offsetX, NaI_Crystal_offsetY, NaI_Crystal_offsetZ);
	G4Transform3D NaI_Crystal_T(RotationNaI, NaI_Crystal_V);

	G4double NaI_CrystalHousingTop_offsetX 	= NaI_Crystal_offsetX;
	G4double NaI_CrystalHousingTop_offsetY 	= NaI_Crystal_offsetY;
	G4double NaI_CrystalHousingTop_offsetZ 	= NaI_Crystal_offsetZ;
	G4ThreeVector NaI_CrystalHousingTop_V (NaI_CrystalHousingTop_offsetX, NaI_CrystalHousingTop_offsetY, NaI_CrystalHousingTop_offsetZ);
	G4Transform3D NaI_CrystalHousingTop_T(RotationNaI, NaI_CrystalHousingTop_V);	

	G4double NaI_CrystalHousingBot_distance = NaI_CrystalHousingTop_height/2 + NaI_CrystalHousingBot_height/2 + 0.1*mm;
	G4double NaI_CrystalHousingBot_offsetX 	= NaI_CrystalHousingTop_offsetX - NaI_CrystalHousingBot_distance*sin(ScatterAngle_rad);
	G4double NaI_CrystalHousingBot_offsetY 	= NaI_CrystalHousingTop_offsetY + NaI_CrystalHousingBot_distance*cos(ScatterAngle_rad);
	G4double NaI_CrystalHousingBot_offsetZ 	= NaI_Crystal_offsetZ;
	G4ThreeVector NaI_CrystalHousingBot_V (NaI_CrystalHousingBot_offsetX, NaI_CrystalHousingBot_offsetY, NaI_CrystalHousingBot_offsetZ);
	G4Transform3D NaI_CrystalHousingBot_T(RotationNaI, NaI_CrystalHousingBot_V);

	G4double NaI_LightShieldBot_distance 	= NaI_CrystalHousingTop_height/2 + NaI_LightShieldBot_height/2;
	G4double NaI_LightShieldBot_offsetX 	= NaI_CrystalHousingTop_offsetX + NaI_LightShieldBot_distance*sin(ScatterAngle_rad);
	G4double NaI_LightShieldBot_offsetY 	= NaI_CrystalHousingTop_offsetY - NaI_LightShieldBot_distance*cos(ScatterAngle_rad);
	G4double NaI_LightShieldBot_offsetZ 	= NaI_Crystal_offsetZ;
	G4ThreeVector NaI_LightShieldBot_V (NaI_LightShieldBot_offsetX, NaI_LightShieldBot_offsetY, NaI_LightShieldBot_offsetZ);
	G4Transform3D NaI_LightShieldBot_T(RotationNaI, NaI_LightShieldBot_V);

	G4double NaI_LightShieldMid_distance 	= NaI_LightShieldBot_height/2 + NaI_LightShieldMid_height/2;
	G4double NaI_LightShieldMid_offsetX 	= NaI_LightShieldBot_offsetX + NaI_LightShieldMid_distance*sin(ScatterAngle_rad);
	G4double NaI_LightShieldMid_offsetY 	= NaI_LightShieldBot_offsetY - NaI_LightShieldMid_distance*cos(ScatterAngle_rad);
	G4double NaI_LightShieldMid_offsetZ 	= NaI_Crystal_offsetZ;
	G4ThreeVector NaI_LightShieldMid_V (NaI_LightShieldMid_offsetX, NaI_LightShieldMid_offsetY, NaI_LightShieldMid_offsetZ);
	G4Transform3D NaI_LightShieldMid_T(RotationNaI, NaI_LightShieldMid_V);

	G4double NaI_LightShieldTop_distance 	= NaI_LightShieldMid_height/2 + NaI_LightShieldTop_height/2;
	G4double NaI_LightShieldTop_offsetX 	= NaI_LightShieldMid_offsetX + NaI_LightShieldTop_distance*sin(ScatterAngle_rad);
	G4double NaI_LightShieldTop_offsetY 	= NaI_LightShieldMid_offsetY - NaI_LightShieldTop_distance*cos(ScatterAngle_rad);
	G4double NaI_LightShieldTop_offsetZ 	= NaI_Crystal_offsetZ;
	G4ThreeVector NaI_LightShieldTop_V (NaI_LightShieldTop_offsetX, NaI_LightShieldTop_offsetY, NaI_LightShieldTop_offsetZ);
	G4Transform3D NaI_LightShieldTop_T(RotationNaI, NaI_LightShieldTop_V);

	// NaI lead castle
	G4double NaIcastle_distance 	= NaI_CrystalHousingBot_height/2 + NaIcastle_offsetDepth - NaIcastle_length/2;
	G4double NaIcastle_offsetX 		= NaI_CrystalHousingBot_offsetX - NaIcastle_distance*sin(ScatterAngle_rad);
	G4double NaIcastle_offsetY 		= NaI_CrystalHousingBot_offsetY + NaIcastle_distance*cos(ScatterAngle_rad);
	G4double NaIcastle_offsetZ 		= NaI_z;
	G4ThreeVector NaIcastle_V (NaIcastle_offsetX, NaIcastle_offsetY, NaIcastle_offsetZ);
	G4Transform3D NaIcastle_T(RotationNaIcastle, NaIcastle_V);

	// lead channel
	G4double LeadChannel_distance 		= NaIcastle_length/2 + LeadChannel_length/2;
	G4double LeadChannel_offsetX 		= NaIcastle_offsetX - LeadChannel_distance*sin(ScatterAngle_rad);
	G4double LeadChannel_offsetY 		= NaIcastle_offsetY + LeadChannel_distance*cos(ScatterAngle_rad);
	G4double LeadChannel_offsetZ 		= NaIcastle_offsetZ;
	G4ThreeVector LeadChannel_V (LeadChannel_offsetX, LeadChannel_offsetY, LeadChannel_offsetZ);
	G4Transform3D LeadChannel_T(RotationNaIcastle, LeadChannel_V);	

	// aperture
	G4double Aperture_distance 	= LeadChannel_length/2 + Aperture_thick/2;
	G4double Aperture_offsetX 	= LeadChannel_offsetX - Aperture_distance*sin(ScatterAngle_rad);
	G4double Aperture_offsetY 	= LeadChannel_offsetY + Aperture_distance*cos(ScatterAngle_rad);
	G4double Aperture_offsetZ 	= LeadChannel_offsetZ;
	G4ThreeVector Aperture_V (Aperture_offsetX, Aperture_offsetY, Aperture_offsetZ);
	G4Transform3D Aperture_T(RotationNaIcastle, Aperture_V);

	// aperture frame
	G4double ApertureFrame_distance = 0.0;
	G4double ApertureFrame_offsetX 	= Aperture_offsetX - ApertureFrame_distance*sin(ScatterAngle_rad);
	G4double ApertureFrame_offsetY 	= Aperture_offsetY + ApertureFrame_distance*cos(ScatterAngle_rad);
	G4double ApertureFrame_offsetZ 	= Aperture_offsetZ;
	G4ThreeVector ApertureFrame_V (ApertureFrame_offsetX, ApertureFrame_offsetY, ApertureFrame_offsetZ);
	G4Transform3D ApertureFrame_T(RotationNaIcastle, ApertureFrame_V);


	// NaI stand, rotated
	G4double NaIstand_x			= NaIcastle_offsetX;
	G4double NaIstand_y			= NaIcastle_offsetY;
	G4double NaIstand_offsetZ	= NaIcastle_height/2;

	G4double NaIstand_PlateTop_offsetX 			= NaIstand_x;
	G4double NaIstand_PlateTop_offsetY 			= NaIstand_y;
	G4double NaIstand_PlateTop_offsetZ 			= NaIcastle_offsetZ - NaIstand_offsetZ - NaIstand_Plate_thick/2;
	G4ThreeVector NaIstand_PlateTop_V (NaIstand_PlateTop_offsetX, NaIstand_PlateTop_offsetY, NaIstand_PlateTop_offsetZ);
	G4Transform3D NaIstand_PlateTop_T(RotationNaIcastle, NaIstand_PlateTop_V);

	G4double NaIstand_VertBar_distance = NaIstand_Plate_side/2-NaIstand_profileOut_out/2;
	
	G4double NaIstand_VertBar_LeftFront_offsetX = NaIstand_x - (NaIstand_VertBar_distance*cos(ScatterAngle_rad)-NaIstand_VertBar_distance*sin(ScatterAngle_rad));
	G4double NaIstand_VertBar_LeftFront_offsetY = NaIstand_y - (NaIstand_VertBar_distance*sin(ScatterAngle_rad)+NaIstand_VertBar_distance*cos(ScatterAngle_rad));
	G4double NaIstand_VertBar_LeftFront_offsetZ = NaIstand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector NaIstand_VertBar_LeftFront_V (NaIstand_VertBar_LeftFront_offsetX, NaIstand_VertBar_LeftFront_offsetY, NaIstand_VertBar_LeftFront_offsetZ);
	G4Transform3D NaIstand_VertBar_LeftFront_T(RotationNaIcastle, NaIstand_VertBar_LeftFront_V);

	G4double NaIstand_VertBar_LeftBack_offsetX 	= NaIstand_x - (NaIstand_VertBar_distance*cos(ScatterAngle_rad)+NaIstand_VertBar_distance*sin(ScatterAngle_rad));
	G4double NaIstand_VertBar_LeftBack_offsetY 	= NaIstand_y - (NaIstand_VertBar_distance*sin(ScatterAngle_rad)-NaIstand_VertBar_distance*cos(ScatterAngle_rad));
	G4double NaIstand_VertBar_LeftBack_offsetZ 	= NaIstand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector NaIstand_VertBar_LeftBack_V (NaIstand_VertBar_LeftBack_offsetX, NaIstand_VertBar_LeftBack_offsetY, NaIstand_VertBar_LeftBack_offsetZ);
	G4Transform3D NaIstand_VertBar_LeftBack_T(RotationNaIcastle, NaIstand_VertBar_LeftBack_V);

	G4double NaIstand_VertBar_RightFront_offsetX = NaIstand_x + (NaIstand_VertBar_distance*cos(ScatterAngle_rad)+NaIstand_VertBar_distance*sin(ScatterAngle_rad));
	G4double NaIstand_VertBar_RightFront_offsetY = NaIstand_y + (NaIstand_VertBar_distance*sin(ScatterAngle_rad)-NaIstand_VertBar_distance*cos(ScatterAngle_rad));
	G4double NaIstand_VertBar_RightFront_offsetZ = NaIstand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector NaIstand_VertBar_RightFront_V (NaIstand_VertBar_RightFront_offsetX, NaIstand_VertBar_RightFront_offsetY, NaIstand_VertBar_RightFront_offsetZ);
	G4Transform3D NaIstand_VertBar_RightFront_T(RotationNaIcastle, NaIstand_VertBar_RightFront_V);

	G4double NaIstand_VertBar_RightBack_offsetX = NaIstand_x + (NaIstand_VertBar_distance*cos(ScatterAngle_rad)-NaIstand_VertBar_distance*sin(ScatterAngle_rad));
	G4double NaIstand_VertBar_RightBack_offsetY = NaIstand_y + (NaIstand_VertBar_distance*sin(ScatterAngle_rad)+NaIstand_VertBar_distance*cos(ScatterAngle_rad));
	G4double NaIstand_VertBar_RightBack_offsetZ = NaIstand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector NaIstand_VertBar_RightBack_V (NaIstand_VertBar_RightBack_offsetX, NaIstand_VertBar_RightBack_offsetY, NaIstand_VertBar_RightBack_offsetZ);
	G4Transform3D NaIstand_VertBar_RightBack_T(RotationNaIcastle, NaIstand_VertBar_RightBack_V);

	G4double NaIstand_PlateBot_offsetX 			= NaIstand_PlateTop_offsetX;
	G4double NaIstand_PlateBot_offsetY 			= NaIstand_PlateTop_offsetY;
	G4double NaIstand_PlateBot_offsetZ 			= NaIstand_VertBar_LeftFront_offsetZ - NaIstand_longProfile_height/2;
	G4ThreeVector NaIstand_PlateBot_V (NaIstand_PlateBot_offsetX, NaIstand_PlateBot_offsetY, NaIstand_PlateBot_offsetZ);
	G4Transform3D NaIstand_PlateBot_T(RotationNaIcastle, NaIstand_PlateBot_V);

	G4double NaIstand_HorzBar_Left_offsetX 		= NaIstand_x - NaIstand_VertBar_distance*cos(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Left_offsetY 		= NaIstand_y - NaIstand_VertBar_distance*sin(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Left_offsetZ 		= NaIstand_VertBar_LeftFront_offsetZ;
	G4ThreeVector NaIstand_HorzBar_Left_V (NaIstand_HorzBar_Left_offsetX, NaIstand_HorzBar_Left_offsetY, NaIstand_HorzBar_Left_offsetZ);
	G4Transform3D NaIstand_HorzBar_Left_T(RotationNaIcastle, NaIstand_HorzBar_Left_V);

	G4double NaIstand_HorzBar_Right_offsetX 	= NaIstand_x + NaIstand_VertBar_distance*cos(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Right_offsetY 	= NaIstand_y + NaIstand_VertBar_distance*sin(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Right_offsetZ 	= NaIstand_VertBar_RightFront_offsetZ;
	G4ThreeVector NaIstand_HorzBar_Right_V (NaIstand_HorzBar_Right_offsetX, NaIstand_HorzBar_Right_offsetY, NaIstand_HorzBar_Right_offsetZ);
	G4Transform3D NaIstand_HorzBar_Right_T(RotationNaIcastle, NaIstand_HorzBar_Right_V);

	G4double NaIstand_HorzBar_Front_offsetX 	= NaIstand_x + NaIstand_VertBar_distance*sin(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Front_offsetY 	= NaIstand_y - NaIstand_VertBar_distance*cos(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Front_offsetZ 	= NaIstand_VertBar_LeftFront_offsetZ;
	G4ThreeVector NaIstand_HorzBar_Front_V (NaIstand_HorzBar_Front_offsetX, NaIstand_HorzBar_Front_offsetY, NaIstand_HorzBar_Front_offsetZ);
	G4Transform3D NaIstand_HorzBar_Front_T(RotationNaIcastle, NaIstand_HorzBar_Front_V);

	G4double NaIstand_HorzBar_Back_offsetX 		= NaIstand_x - NaIstand_VertBar_distance*sin(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Back_offsetY 		= NaIstand_y + NaIstand_VertBar_distance*cos(ScatterAngle_rad);
	G4double NaIstand_HorzBar_Back_offsetZ 		= NaIstand_VertBar_LeftBack_offsetZ;
	G4ThreeVector NaIstand_HorzBar_Back_V (NaIstand_HorzBar_Back_offsetX, NaIstand_HorzBar_Back_offsetY, NaIstand_HorzBar_Back_offsetZ);
	G4Transform3D NaIstand_HorzBar_Back_T(RotationNaIcastle, NaIstand_HorzBar_Back_V);

	G4double NaIstand_CradleBar_distance 		= NaIstand_CradlePanel_length/2 - NaIstand_CradleBar_oD/2;

	G4double NaIstand_CradleBarLeft_offsetX 	= NaIstand_x - NaIstand_CradleBar_distance*cos(ScatterAngle_rad);
	G4double NaIstand_CradleBarLeft_offsetY 	= NaIstand_y - NaIstand_CradleBar_distance*sin(ScatterAngle_rad);
	G4double NaIstand_CradleBarLeft_offsetZ 	= NaIstand_PlateTop_offsetZ + NaIstand_Plate_thick/2 + NaIstand_CradleBar_height/2;
	G4ThreeVector NaIstand_CradleBarLeft_V (NaIstand_CradleBarLeft_offsetX, NaIstand_CradleBarLeft_offsetY, NaIstand_CradleBarLeft_offsetZ);
	G4Transform3D NaIstand_CradleBarLeft_T(RotationNaIcastle, NaIstand_CradleBarLeft_V);

	G4double NaIstand_CradleBarRight_offsetX 	= NaIstand_x + NaIstand_CradleBar_distance*cos(ScatterAngle_rad);
	G4double NaIstand_CradleBarRight_offsetY 	= NaIstand_y + NaIstand_CradleBar_distance*sin(ScatterAngle_rad);
	G4double NaIstand_CradleBarRight_offsetZ 	= NaIstand_PlateTop_offsetZ + NaIstand_Plate_thick/2 + NaIstand_CradleBar_height/2;
	G4ThreeVector NaIstand_CradleBarRight_V (NaIstand_CradleBarRight_offsetX, NaIstand_CradleBarRight_offsetY, NaIstand_CradleBarRight_offsetZ);
	G4Transform3D NaIstand_CradleBarRight_T(RotationNaIcastle, NaIstand_CradleBarRight_V);

	G4double NaIstand_CradlePanel_offsetX 		= NaIstand_x;
	G4double NaIstand_CradlePanel_offsetY 		= NaIstand_y;
	G4double NaIstand_CradlePanel_offsetZ 		= NaIstand_CradleBarLeft_offsetZ + NaIstand_CradleBar_height/2 + NaIstand_CradlePanel_thick/2;
	G4ThreeVector NaIstand_CradlePanel_V (NaIstand_CradlePanel_offsetX, NaIstand_CradlePanel_offsetY, NaIstand_CradlePanel_offsetZ);
	G4Transform3D NaIstand_CradlePanel_T(RotationNaIcastle, NaIstand_CradlePanel_V);

	G4cout <<"----------------------------------------------"<< G4endl;
	G4cout <<"| > NaI distance:      "<< NaI_distance 		<<" mm"<< G4endl;
	G4cout <<"| > Scattering angle:  "<< ScatterAngle_deg 	<<" deg  =  "<< ScatterAngle_rad <<" rad"<< G4endl;
	G4cout <<"| > X:                 "<< NaI_x 				<<" mm"<< G4endl;
	G4cout <<"| > Y:                 "<< NaI_y 				<<" mm"<< G4endl;
	G4cout <<"| > Z:                 "<< NaI_z 				<<" mm"<< G4endl;
	G4cout <<"----------------------------------------------"<< G4endl;

	// Lead Channel stand, rotated
	G4double LeadChannelStand_x			= LeadChannel_offsetX;
	G4double LeadChannelStand_y			= LeadChannel_offsetY;
	G4double LeadChannelStand_offsetZ	= LeadChannel_height/2;

	G4double LeadChannelStand_PlateTop_offsetX 			= LeadChannelStand_x;
	G4double LeadChannelStand_PlateTop_offsetY 			= LeadChannelStand_y;
	G4double LeadChannelStand_PlateTop_offsetZ 			= LeadChannel_offsetZ - LeadChannelStand_offsetZ - NaIstand_Plate_thick/2;
	G4ThreeVector LeadChannelStand_PlateTop_V (LeadChannelStand_PlateTop_offsetX, LeadChannelStand_PlateTop_offsetY, LeadChannelStand_PlateTop_offsetZ);
	G4Transform3D LeadChannelStand_PlateTop_T(RotationNaIcastle, LeadChannelStand_PlateTop_V);

	G4double LeadChannelStand_VertBar_distance = NaIstand_Plate_side/2-NaIstand_profileOut_out/2;
	
	G4double LeadChannelStand_VertBar_LeftFront_offsetX = LeadChannelStand_x - (LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad)-LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_LeftFront_offsetY = LeadChannelStand_y - (LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad)+LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_LeftFront_offsetZ = LeadChannelStand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector LeadChannelStand_VertBar_LeftFront_V (LeadChannelStand_VertBar_LeftFront_offsetX, LeadChannelStand_VertBar_LeftFront_offsetY, LeadChannelStand_VertBar_LeftFront_offsetZ);
	G4Transform3D LeadChannelStand_VertBar_LeftFront_T(RotationNaIcastle, LeadChannelStand_VertBar_LeftFront_V);

	G4double LeadChannelStand_VertBar_LeftBack_offsetX 	= LeadChannelStand_x - (LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad)+LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_LeftBack_offsetY 	= LeadChannelStand_y - (LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad)-LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_LeftBack_offsetZ 	= LeadChannelStand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector LeadChannelStand_VertBar_LeftBack_V (LeadChannelStand_VertBar_LeftBack_offsetX, LeadChannelStand_VertBar_LeftBack_offsetY, LeadChannelStand_VertBar_LeftBack_offsetZ);
	G4Transform3D LeadChannelStand_VertBar_LeftBack_T(RotationNaIcastle, LeadChannelStand_VertBar_LeftBack_V);

	G4double LeadChannelStand_VertBar_RightFront_offsetX = LeadChannelStand_x + (LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad)+LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_RightFront_offsetY = LeadChannelStand_y + (LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad)-LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_RightFront_offsetZ = LeadChannelStand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector LeadChannelStand_VertBar_RightFront_V (LeadChannelStand_VertBar_RightFront_offsetX, LeadChannelStand_VertBar_RightFront_offsetY, LeadChannelStand_VertBar_RightFront_offsetZ);
	G4Transform3D LeadChannelStand_VertBar_RightFront_T(RotationNaIcastle, LeadChannelStand_VertBar_RightFront_V);

	G4double LeadChannelStand_VertBar_RightBack_offsetX = LeadChannelStand_x + (LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad)-LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_RightBack_offsetY = LeadChannelStand_y + (LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad)+LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad));
	G4double LeadChannelStand_VertBar_RightBack_offsetZ = LeadChannelStand_PlateTop_offsetZ - NaIstand_Plate_thick/2 - NaIstand_longProfile_height/2;
	G4ThreeVector LeadChannelStand_VertBar_RightBack_V (LeadChannelStand_VertBar_RightBack_offsetX, LeadChannelStand_VertBar_RightBack_offsetY, LeadChannelStand_VertBar_RightBack_offsetZ);
	G4Transform3D LeadChannelStand_VertBar_RightBack_T(RotationNaIcastle, LeadChannelStand_VertBar_RightBack_V);

	G4double LeadChannelStand_PlateBot_offsetX 			= LeadChannelStand_PlateTop_offsetX;
	G4double LeadChannelStand_PlateBot_offsetY 			= LeadChannelStand_PlateTop_offsetY;
	G4double LeadChannelStand_PlateBot_offsetZ 			= LeadChannelStand_VertBar_LeftFront_offsetZ - NaIstand_longProfile_height/2;
	G4ThreeVector LeadChannelStand_PlateBot_V (LeadChannelStand_PlateBot_offsetX, LeadChannelStand_PlateBot_offsetY, LeadChannelStand_PlateBot_offsetZ);
	G4Transform3D LeadChannelStand_PlateBot_T(RotationNaIcastle, LeadChannelStand_PlateBot_V);

	G4double LeadChannelStand_HorzBar_Left_offsetX 		= LeadChannelStand_x - LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Left_offsetY 		= LeadChannelStand_y - LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Left_offsetZ 		= LeadChannelStand_VertBar_LeftFront_offsetZ;
	G4ThreeVector LeadChannelStand_HorzBar_Left_V (LeadChannelStand_HorzBar_Left_offsetX, LeadChannelStand_HorzBar_Left_offsetY, LeadChannelStand_HorzBar_Left_offsetZ);
	G4Transform3D LeadChannelStand_HorzBar_Left_T(RotationNaIcastle, LeadChannelStand_HorzBar_Left_V);

	G4double LeadChannelStand_HorzBar_Right_offsetX 	= LeadChannelStand_x + LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Right_offsetY 	= LeadChannelStand_y + LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Right_offsetZ 	= LeadChannelStand_VertBar_RightFront_offsetZ;
	G4ThreeVector LeadChannelStand_HorzBar_Right_V (LeadChannelStand_HorzBar_Right_offsetX, LeadChannelStand_HorzBar_Right_offsetY, LeadChannelStand_HorzBar_Right_offsetZ);
	G4Transform3D LeadChannelStand_HorzBar_Right_T(RotationNaIcastle, LeadChannelStand_HorzBar_Right_V);

	G4double LeadChannelStand_HorzBar_Front_offsetX 	= LeadChannelStand_x + LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Front_offsetY 	= LeadChannelStand_y - LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Front_offsetZ 	= LeadChannelStand_VertBar_LeftFront_offsetZ;
	G4ThreeVector LeadChannelStand_HorzBar_Front_V (LeadChannelStand_HorzBar_Front_offsetX, LeadChannelStand_HorzBar_Front_offsetY, LeadChannelStand_HorzBar_Front_offsetZ);
	G4Transform3D LeadChannelStand_HorzBar_Front_T(RotationNaIcastle, LeadChannelStand_HorzBar_Front_V);

	G4double LeadChannelStand_HorzBar_Back_offsetX 		= LeadChannelStand_x - LeadChannelStand_VertBar_distance*sin(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Back_offsetY 		= LeadChannelStand_y + LeadChannelStand_VertBar_distance*cos(ScatterAngle_rad);
	G4double LeadChannelStand_HorzBar_Back_offsetZ 		= LeadChannelStand_VertBar_LeftBack_offsetZ;
	G4ThreeVector LeadChannelStand_HorzBar_Back_V (LeadChannelStand_HorzBar_Back_offsetX, LeadChannelStand_HorzBar_Back_offsetY, LeadChannelStand_HorzBar_Back_offsetZ);
	G4Transform3D LeadChannelStand_HorzBar_Back_T(RotationNaIcastle, LeadChannelStand_HorzBar_Back_V);

	// collimator
	//G4double Collimator_distance 	= 70 *cm;
	G4double Collimator_distance 	= 33 *cm; // 28 cm 'LXe Center - collimator exit'
	G4double Collimator_x 			= Xuerich_x;
	G4double Collimator_y 			= Xuerich_y + Collimator_distance;
	G4double Collimator_z 			= NaI_z;
	G4double CollimatorOpening_y 	= Collimator_y - CollimatorBox_height/2;
	G4ThreeVector Collimator_V (Collimator_x, Collimator_y, Collimator_z);
	G4Transform3D Collimator_T(RotationXPlus90, Collimator_V);
	G4cout <<"| > Collimator distance: "<< Collimator_distance 	<<" mm"<< G4endl;
	G4cout <<"| > X:                   "<< Collimator_x 		<<" mm"<< G4endl;
	G4cout <<"| > Y:                   "<< Collimator_y 		<<" mm"<< G4endl;
	G4cout <<"| > Z:                   "<< Collimator_z 		<<" mm"<< G4endl;
	G4cout <<"| > Opening at:          "<< CollimatorOpening_y 	<<" mm"<< G4endl;
	G4cout <<"----------------------------------------------"<< G4endl;


	//==============================================================================================
	//==============	Physical Volumes (declared in 'XuerichDetectorGeometry.hh')	================
	//==============================================================================================
	G4PVPlacement *Laboratory_phys					= new G4PVPlacement(0, G4ThreeVector(0, 0, 0), Laboratory_log, "Laboratory_phys", 0, false, 0);

	// XuStand, ROTATED
	G4PVPlacement *XuStand_HorzBar_TopLeft_phys		= new G4PVPlacement(XuStand_HorzBar_TopLeft_T, XuStand_HorzBar_TopLeft_log, "XuStand_HorzBar_TopLeft_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_TopRight_phys	= new G4PVPlacement(XuStand_HorzBar_TopRight_T, XuStand_HorzBar_TopRight_log, "XuStand_HorzBar_TopRight_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_TopFront_phys	= new G4PVPlacement(XuStand_HorzBar_TopFront_T, XuStand_HorzBar_TopFront_log, "XuStand_HorzBar_TopFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_TopBack_phys		= new G4PVPlacement(XuStand_HorzBar_TopBack_T, XuStand_HorzBar_TopBack_log, "XuStand_HorzBar_TopBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_WingFront_phys			= new G4PVPlacement(XuStand_WingFront_T, XuStand_WingFront_log, "XuStand_WingFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_WingBack_phys			= new G4PVPlacement(XuStand_WingBack_T, XuStand_WingBack_log, "XuStand_WingBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_MidLeft_phys		= new G4PVPlacement(XuStand_HorzBar_MidLeft_T, XuStand_HorzBar_MidLeft_log, "XuStand_HorzBar_MidLeft_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_MidRight_phys	= new G4PVPlacement(XuStand_HorzBar_MidRight_T, XuStand_HorzBar_MidRight_log, "XuStand_HorzBar_MidRight_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_MidFront_phys	= new G4PVPlacement(XuStand_HorzBar_MidFront_T, XuStand_HorzBar_MidFront_log, "XuStand_HorzBar_MidFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_MidBack_phys		= new G4PVPlacement(XuStand_HorzBar_MidBack_T, XuStand_HorzBar_MidBack_log, "XuStand_HorzBar_MidBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_VertBar_LeftFront_phys	= new G4PVPlacement(XuStand_VertBar_LeftFront_T, XuStand_VertBar_LeftFront_log, "XuStand_VertBar_LeftFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_VertBar_LeftBack_phys	= new G4PVPlacement(XuStand_VertBar_LeftBack_T, XuStand_VertBar_LeftBack_log, "XuStand_VertBar_LeftBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_VertBar_RightFront_phys	= new G4PVPlacement(XuStand_VertBar_RightFront_T, XuStand_VertBar_RightFront_log, "XuStand_VertBar_RightFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_VertBar_RightBack_phys	= new G4PVPlacement(XuStand_VertBar_RightBack_T, XuStand_VertBar_RightBack_log, "XuStand_VertBar_RightBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_BotLeft_phys		= new G4PVPlacement(XuStand_HorzBar_BotLeft_T, XuStand_HorzBar_BotLeft_log, "XuStand_HorzBar_BotLeft_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_BotRight_phys	= new G4PVPlacement(XuStand_HorzBar_BotRight_T, XuStand_HorzBar_BotRight_log, "XuStand_HorzBar_BotRight_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_BotFront_phys	= new G4PVPlacement(XuStand_HorzBar_BotFront_T, XuStand_HorzBar_BotFront_log, "XuStand_HorzBar_BotFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_HorzBar_BotBack_phys		= new G4PVPlacement(XuStand_HorzBar_BotBack_T, XuStand_HorzBar_BotBack_log, "XuStand_HorzBar_BotBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_BottomPlate_phys			= new G4PVPlacement(XuStand_BottomPlate_T, XuStand_BottomPlate_log, "XuStand_BottomPlate_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CylinderLeftFront_phys	= new G4PVPlacement(XuStand_CylinderLeftFront_T, XuStand_CylinderLeftFront_log, "XuStand_CylinderLeftFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CylinderLeftBack_phys	= new G4PVPlacement(XuStand_CylinderLeftBack_T, XuStand_CylinderLeftBack_log, "XuStand_CylinderLeftBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CylinderRightFront_phys	= new G4PVPlacement(XuStand_CylinderRightFront_T, XuStand_CylinderRightFront_log, "XuStand_CylinderRightFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CylinderRightBack_phys	= new G4PVPlacement(XuStand_CylinderRightBack_T, XuStand_CylinderRightBack_log, "XuStand_CylinderRightBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CradleBarLeftFront_phys	= new G4PVPlacement(XuStand_CradleBarLeftFront_T, XuStand_CradleBarLeftFront_log, "XuStand_CradleBarLeftFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CradleBarLeftBack_phys	= new G4PVPlacement(XuStand_CradleBarLeftBack_T, XuStand_CradleBarLeftBack_log, "XuStand_CradleBarLeftBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CradleBarRightFront_phys	= new G4PVPlacement(XuStand_CradleBarRightFront_T, XuStand_CradleBarRightFront_log, "XuStand_CradleBarRightFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CradleBarRightBack_phys	= new G4PVPlacement(XuStand_CradleBarRightBack_T, XuStand_CradleBarRightBack_log, "XuStand_CradleBarRightBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CradlePanelLong_phys		= new G4PVPlacement(XuStand_CradlePanelLong_T, XuStand_CradlePanelLong_log, "XuStand_CradlePanelLong_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_CradlePanelShort_phys	= new G4PVPlacement(XuStand_CradlePanelShort_T, XuStand_CradlePanelShort_log, "XuStand_CradlePanelShort_phys", Laboratory_log, false, 0);
	G4PVPlacement *XuStand_Goniometer_phys			= new G4PVPlacement(0, G4ThreeVector(XuStand_Goniometer_offsetX, XuStand_Goniometer_offsetY, XuStand_Goniometer_offsetZ), XuStand_Goniometer_log, "XuStand_Goniometer_phys", Laboratory_log, false, 0);

	// dewar
	G4PVPlacement *DewarOutTube_phys				= new G4PVPlacement(0, G4ThreeVector(DewarOutTube_offsetX, DewarOutTube_offsetY, DewarOutTube_offsetZ), DewarOutTube_log, "DewarOutTube_phys", Laboratory_log, false, 0);
	G4PVPlacement *DewarOutBot_phys					= new G4PVPlacement(0, G4ThreeVector(DewarOutBot_offsetX, DewarOutBot_offsetY, DewarOutBot_offsetZ), DewarOutBot_log, "DewarOutBot_phys", Laboratory_log, false, 0);
	G4PVPlacement *DewarInTube_phys					= new G4PVPlacement(0, G4ThreeVector(DewarInTube_offsetX, DewarInTube_offsetY, DewarInTube_offsetZ), DewarInTube_log, "DewarInTube_phys", Laboratory_log, false, 0);
	G4PVPlacement *DewarInBot_phys					= new G4PVPlacement(0, G4ThreeVector(DewarInBot_offsetX, DewarInBot_offsetY, DewarInBot_offsetZ), DewarInBot_log, "DewarInBot_phys", Laboratory_log, false, 0);
	G4PVPlacement *DewarTop_phys					= new G4PVPlacement(0, G4ThreeVector(DewarTop_offsetX, DewarTop_offsetY, DewarTop_offsetZ), DewarTop_log, "DewarTop_phys", Laboratory_log, false, 0);
	// outer cryostat can
	G4PVPlacement *OuterCanUpperFlange_phys			= new G4PVPlacement(0, G4ThreeVector(Xuerich_x, Xuerich_y, OuterCanUpperFlange_offsetZ), OuterCanUpperFlange_log, "OuterCanUpperFlange_phys", Laboratory_log, false, 0);
	G4PVPlacement *OuterCanLowerFlange_phys			= new G4PVPlacement(0, G4ThreeVector(Xuerich_x, Xuerich_y, OuterCanLowerFlange_offsetZ), OuterCanLowerFlange_log, "OuterCanLowerFlange_phys", Laboratory_log, false, 0);
	G4PVPlacement *OuterCanTube_phys				= new G4PVPlacement(0, G4ThreeVector(Xuerich_x, Xuerich_y, OuterCanTube_offsetZ), OuterCanTube_log, "OuterCanTube_phys", Laboratory_log, false, 0);
	G4PVPlacement *OuterCanBottom_phys				= new G4PVPlacement(0, G4ThreeVector(Xuerich_x, Xuerich_y, OuterCanBottom_offsetZ),	OuterCanBottom_log, "OuterCanBottom_phys", Laboratory_log, false, 0);
	// vacuum in the cryostat
	G4PVPlacement *CryostatVacuum_phys				= new G4PVPlacement(0, G4ThreeVector(Xuerich_x, Xuerich_y, CryostatVacuum_offsetZ), CryostatVacuum_log, "CryostatVacuum_phys",	Laboratory_log,	false, 0);
	// feedthroughs inside the outer cryostat
	G4PVPlacement *Pipe1_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe1_offsetX, Pipe1_offsetY, PipeIn_offsetZ), Pipe1in_log, "Pipe1in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe2_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe2_offsetX, Pipe2_offsetY, PipeIn_offsetZ), Pipe2in_log, "Pipe2in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe3_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe3_offsetX, Pipe3_offsetY, PipeIn_offsetZ), Pipe3in_log, "Pipe3in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe4_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe4_offsetX, Pipe4_offsetY, PipeIn_offsetZ), Pipe4in_log, "Pipe4in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe5_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe5_offsetX, Pipe5_offsetY, PipeIn_offsetZ), Pipe5in_log, "Pipe5in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe6_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe6_offsetX, Pipe6_offsetY, PipeIn_offsetZ), Pipe6in_log, "Pipe6in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe7_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe7_offsetX, Pipe7_offsetY, PipeIn_offsetZ), Pipe7in_log, "Pipe7in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe8_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe8_offsetX, Pipe8_offsetY, PipeIn_offsetZ), Pipe8in_log, "Pipe8in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe9_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe9_offsetX, Pipe9_offsetY, PipeIn_offsetZ), Pipe9in_log, "Pipe9in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe10_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe10_offsetX, Pipe10_offsetY, PipeIn_offsetZ), Pipe10in_log, "Pipe10in_phys", CryostatVacuum_log, false, 0);	
	G4PVPlacement *Pipe11_phys						= new G4PVPlacement(0, G4ThreeVector(Pipe11_offsetX, Pipe11_offsetY, PipeIn_offsetZ), Pipe11in_log, "Pipe11in_phys", CryostatVacuum_log, false, 0);	
	// feedthroughs outside
	G4PVPlacement *Pipe1out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe1_offsetX, Pipe1_offsetY, Pipe1out_offsetZ), Pipe1out_log, "Pipe1out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe2out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe2_offsetX, Pipe2_offsetY, Pipe2out_offsetZ), Pipe2out_log, "Pipe2out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe3out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe3_offsetX, Pipe3_offsetY, Pipe3out_offsetZ), Pipe3out_log, "Pipe3out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe4out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe4_offsetX, Pipe4_offsetY, Pipe4out_offsetZ), Pipe4out_log, "Pipe4out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe5out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe5_offsetX, Pipe5_offsetY, Pipe5out_offsetZ), Pipe5out_log, "Pipe5out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe6out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe6_offsetX, Pipe6_offsetY, Pipe6out_offsetZ), Pipe6out_log, "Pipe6out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe7out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe7_offsetX, Pipe7_offsetY, Pipe7out_offsetZ), Pipe7out_log, "Pipe7out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe8out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe8_offsetX, Pipe8_offsetY, Pipe8out_offsetZ), Pipe8out_log, "Pipe8out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe9out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe9_offsetX, Pipe9_offsetY, Pipe9out_offsetZ), Pipe9out_log, "Pipe9out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe10out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe10_offsetX, Pipe10_offsetY, Pipe10out_offsetZ), Pipe10out_log, "Pipe10out_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe11out_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe11_offsetX, Pipe11_offsetY, Pipe11out_offsetZ), Pipe11out_log, "Pipe11out_phys", Laboratory_log, false, 0);	

	G4PVPlacement *Pipe1flange_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe1_offsetX, Pipe1_offsetY, Pipe1flange_offsetZ), Pipe1flange_log, "Pipe1flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe2flange_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe2_offsetX, Pipe2_offsetY, Pipe2flange_offsetZ), Pipe2flange_log, "Pipe2flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe3flange_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe3_offsetX, Pipe3_offsetY, Pipe3flange_offsetZ), Pipe3flange_log, "Pipe3flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe4flange_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe4_offsetX, Pipe4_offsetY, Pipe4flange_offsetZ), Pipe4flange_log, "Pipe4flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe5flange1_phys				= new G4PVPlacement(0, G4ThreeVector(Pipe5_offsetX, Pipe5_offsetY, Pipe5flange1_offsetZ), Pipe5flange1_log, "Pipe5flange1_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe5flange2_phys				= new G4PVPlacement(0, G4ThreeVector(Pipe5_offsetX, Pipe5_offsetY, Pipe5flange2_offsetZ), Pipe5flange1_log, "Pipe5flange2_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe6flange_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe6_offsetX, Pipe6_offsetY, Pipe6flange_offsetZ), Pipe6flange_log, "Pipe6flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe8flange_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe8_offsetX, Pipe8_offsetY, Pipe8flange_offsetZ), Pipe8flange_log, "Pipe8flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe9flange_phys					= new G4PVPlacement(0, G4ThreeVector(Pipe9_offsetX, Pipe9_offsetY, Pipe9flange_offsetZ), Pipe9flange_log, "Pipe9flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe10flange_phys				= new G4PVPlacement(0, G4ThreeVector(Pipe10_offsetX, Pipe10_offsetY, Pipe10flange_offsetZ), Pipe10flange_log, "Pipe10flange_phys", Laboratory_log, false, 0);	
	G4PVPlacement *Pipe11flange_phys				= new G4PVPlacement(0, G4ThreeVector(Pipe11_offsetX, Pipe11_offsetY, Pipe11flange_offsetZ), Pipe11flange_log, "Pipe11flange_phys", Laboratory_log, false, 0);	
	// inner cryostat can
	G4PVPlacement *InnerCanTube_phys				= new G4PVPlacement(0, G4ThreeVector(0, 0, InnerCanTube_offsetZ), InnerCanTube_log, "InnerCanTube_phys", CryostatVacuum_log, false, 0);
	G4PVPlacement *InnerCanBottom_phys				= new G4PVPlacement(0, G4ThreeVector(0, 0, InnerCanBottom_offsetZ),	InnerCanBottom_log, "InnerCanBottom_phys", CryostatVacuum_log, false, 0);
	G4PVPlacement *InnerCanLowerFlange_phys			= new G4PVPlacement(0, G4ThreeVector(0, 0, InnerCanLowerFlange_offsetZ), InnerCanLowerFlange_log, "InnerCanLowerFlange_phys", CryostatVacuum_log, false, 0);
	G4PVPlacement *InnerCanUpperFlange_phys			= new G4PVPlacement(0, G4ThreeVector(0, 0, InnerCanUpperFlange_offsetZ), InnerCanUpperFlange_log, "InnerCanUpperFlange_phys", CryostatVacuum_log, false, 0);
	// Aluminum can
	G4PVPlacement *AlCanTube_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AlCanTube_offsetZ), AlCanTube_log, "AlCanTube_phys", CryostatVacuum_log, false, 0);
	G4PVPlacement *AlCanBottom_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AlCanBottom_offsetZ), AlCanBottom_log, "AlCanBottom_phys", CryostatVacuum_log, false, 0);
	// copper finger
	G4PVPlacement *CopperFinger_phys				= new G4PVPlacement(0, G4ThreeVector(CopperFinger_offsetX, CopperFinger_offsetY, CopperFinger_offsetZ), CopperFinger_log, "CopperFinger_phys", Laboratory_log, false, 0);
	// xenon
	G4PVPlacement *LXeVol_phys						= new G4PVPlacement(0, G4ThreeVector(0, 0, LXeVol_offsetZ),	LXeVol_log, "LXeVol_phys", CryostatVacuum_log, false, 0);
	G4PVPlacement *GXeVol_phys						= new G4PVPlacement(0, G4ThreeVector(0, 0, GXeVol_offsetZ),	GXeVol_log, "GXeVol_phys", CryostatVacuum_log, false, 0);
	G4PVPlacement *LXeTarget_phys					= new G4PVPlacement(0, G4ThreeVector(0, 0, LXeTarget_offsetZ),	LXeTarget_log, "LXeTarget_phys", LXeVol_log, false, 0);
	G4PVPlacement *LXeGate_phys						= new G4PVPlacement(0, G4ThreeVector(0, 0, LXeGate_offsetZ),	LXeGate_log, "LXeGate_phys", GXeVol_log, false, 0);
	//G4PVPlacement *LXeCathode_phys				= new G4PVPlacement(0, G4ThreeVector(0, 0, LXeCathode_offsetZ),	LXeCathode_log, "LXeCathode_phys", LXeVol_log, false, 0);
	// parts in GXe
	G4PVPlacement *Washer1_phys						= new G4PVPlacement(0, G4ThreeVector(Washer1_offsetX, Washer1_offsetY, Washer_offsetZ), Washer1_log, "Washer1_phys", GXeVol_log, false, 0);
	G4PVPlacement *Washer2_phys						= new G4PVPlacement(0, G4ThreeVector(Washer2_offsetX, Washer2_offsetY, Washer_offsetZ), Washer2_log, "Washer2_phys", GXeVol_log, false, 0);
	G4PVPlacement *Washer3_phys						= new G4PVPlacement(0, G4ThreeVector(Washer3_offsetX, Washer3_offsetY, Washer_offsetZ), Washer3_log, "Washer3_phys", GXeVol_log, false, 0);
	G4PVPlacement *Washer4_phys						= new G4PVPlacement(0, G4ThreeVector(Washer4_offsetX, Washer4_offsetY, Washer_offsetZ), Washer4_log, "Washer4_phys", GXeVol_log, false, 0);
	G4PVPlacement *SteelHolderUpperSquare_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, SteelHolderUpperSquare_offsetZ), SteelHolderUpperSquare_log,	"SteelHolderUpperSquare_phys", GXeVol_log, false, 0);
	G4PVPlacement *SteelHolder_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, SteelHolderTube_offsetZ), SteelHolderTube_log, "SteelHolderTube_phys", GXeVol_log, false, 0);	
	G4PVPlacement *SteelHolderLowerSquare_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, SteelHolderLowerSquare_offsetZ),	SteelHolderLowerSquare_log,	"SteelHolderLowerSquare_phys", GXeVol_log, false, 0);
	G4PVPlacement *Support1_phys					= new G4PVPlacement(0, G4ThreeVector(Support1_offsetX, Support1_offsetY, Support_offsetZ), Support1_log, "Support1_phys", GXeVol_log, false, 0);
	G4PVPlacement *Support2_phys					= new G4PVPlacement(0, G4ThreeVector(Support2_offsetX, Support2_offsetY, Support_offsetZ), Support2_log, "Support2_phys", GXeVol_log, false, 0);
	G4PVPlacement *Support3_phys					= new G4PVPlacement(0, G4ThreeVector(Support3_offsetX, Support3_offsetY, Support_offsetZ), Support1_log, "Support3_phys", GXeVol_log, false, 0);
	G4PVPlacement *Support4_phys					= new G4PVPlacement(0, G4ThreeVector(Support4_offsetX,	Support4_offsetY,	Support_offsetZ),	Support4_log, "Support4_phys", GXeVol_log, false, 0);
	G4PVPlacement *TopClamp_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, TopClamp_offsetZ), TopClamp_log, "TopClamp_phys", GXeVol_log, false, 0);
	G4PVPlacement *TopHolderUpperTube_phys			= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, TopHolderUpperTube_offsetZ),	TopHolderUpperTube_log, "TopHolderUpperTube_phys", GXeVol_log, false, 0);
	G4PVPlacement *TopHolderLowerTube_phys			= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, TopHolderLowerTube_offsetZ),	TopHolderLowerTube_log, "TopHolderLowerTube_phys", GXeVol_log, false, 0);
	G4PVPlacement *TopHolderSquare_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, TopHolderSquare_offsetZ), TopHolderSquare_log, "TopHolderSquare_phys", GXeVol_log, false, 0);
	G4PVPlacement *AnodeGrid_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, AnodeGrid_offsetZ), AnodeGrid_log, "AnodeGrid_phys", GXeVol_log, false, 0);
	G4PVPlacement *ExtractionSpacerUpperHalf_phys	= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, ExtractionSpacerUpperHalf_offsetZ), ExtractionSpacerUpperHalf_log, "ExtractionSpacerUpperHalf_phys",	GXeVol_log, false, 0);
	// parts in LXe
	G4PVPlacement *ExtractionSpacerLowerHalf_phys	= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, ExtractionSpacerLowerHalf_offsetZ), ExtractionSpacerLowerHalf_log,	"ExtractionSpacerLowerHalf_phys",	LXeVol_log,		false, 0);
	G4PVPlacement *GateGrid_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, GateGrid_offsetZ), GateGrid_log, "GateGrid_phys", LXeVol_log, false, 0);
	G4PVPlacement *DriftSpacerUpperSquare_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, DriftSpacerUpperSquare_offsetZ),	DriftSpacerUpperSquare_log, "DriftSpacerUpperSquare_phys", LXeVol_log, false, 0);
	G4PVPlacement *DriftSpacerUpperTube_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, DriftSpacerUpperTube_offsetZ), DriftSpacerUpperTube_log, "DriftSpacerUpperTube_phys", LXeVol_log, false, 0);
	G4PVPlacement *DriftSpacerMiddleTube_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, DriftSpacerMiddleTube_offsetZ), DriftSpacerMiddleTube_log, "DriftSpacerMiddleTube_phys", LXeVol_log, false, 0);	
	G4PVPlacement *DriftSpacerLowerTube_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, DriftSpacerLowerTube_offsetZ), DriftSpacerLowerTube_log, "DriftSpacerLowerTube_phys", LXeVol_log, false, 0);
	G4PVPlacement *DriftSpacerLowerSquare_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, DriftSpacerLowerSquare_offsetZ), DriftSpacerLowerSquare_log, "DriftSpacerLowerSquare_phys", LXeVol_log, false, 0);
	G4PVPlacement *CathodeGrid_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, CathodeGrid_offsetZ), CathodeGrid_log, "CathodeGrid_phys", LXeVol_log, false, 0);
	G4PVPlacement *BottomHolderUpperTube_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, BottomHolderUpperTube_offsetZ), BottomHolderUpperTube_log, "BottomHolderUpperTube_phys", LXeVol_log, false, 0);
	G4PVPlacement *BottomHolderMiddleTube_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, BottomHolderMiddleTube_offsetZ), BottomHolderMiddleTube_log, "BottomHolderMiddleTube_phys", LXeVol_log, false, 0);
	G4PVPlacement *BottomHolderLowerTube_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, BottomHolderLowerTube_offsetZ), BottomHolderLowerTube_log, "BottomHolderLowerTube_phys", LXeVol_log, false, 0);
	G4PVPlacement *Filler_phys						= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, Filler_offsetZ),	Filler_log, "Filler_phys", LXeVol_log, false, 0);
	G4PVPlacement *BottomClampUpperPart_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, BottomClampUpperPart_offsetZ), BottomClampUpperPart_log, "BottomClampUpperPart_phys", LXeVol_log, false, 0);
	G4PVPlacement *BottomClampMiddlePart_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, BottomClampMiddlePart_offsetZ), BottomClampMiddlePart_log, "BottomClampMiddlePart_phys", LXeVol_log, false, 0);
	G4PVPlacement *BottomClampLowerPart_phys		= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, BottomClampLowerPart_offsetZ), BottomClampLowerPart_log, "BottomClampLowerPart_phys", LXeVol_log, false, 0);
	// PMTs
	G4PVPlacement *topPMTcasing_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, topPMTcasing_offsetZ), topPMTcasing_log, "topPMTcasing_phys", GXeVol_log, false, 0);
	G4PVPlacement *topPMTinterior_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, topPMTinterior_offsetZ), topPMTinterior_log, "topPMTinterior_phys", topPMTcasing_log, false, 0);
	G4PVPlacement *topPMTwindow_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, topPMTwindow_offsetZ), topPMTwindow_log, "topPMTwindow_phys", topPMTinterior_log, false, 0);
	G4PVPlacement *topPMTcathode_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, topPMTcathode_offsetZ), topPMTcathode_log, "topPMTcathode_phys",	topPMTwindow_log, false, 0);
	G4PVPlacement *bottomPMTcasing_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bottomPMTcasing_offsetZ), bottomPMTcasing_log, "bottomPMTcasing_phys", LXeVol_log, false, 0);
	G4PVPlacement *bottomPMTinterior_phys			= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bottomPMTinterior_offsetZ), bottomPMTinterior_log, "bottomPMTinterior_phys",	bottomPMTcasing_log, false, 0);
	G4PVPlacement *bottomPMTwindow_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bottomPMTwindow_offsetZ), bottomPMTwindow_log, "bottomPMTwindow_phys", bottomPMTinterior_log, false, 0);
	G4PVPlacement *bottomPMTcathode_phys			= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bottomPMTcathode_offsetZ), bottomPMTcathode_log, "bottomPMTcathode_phys", bottomPMTwindow_log, false, 0);
	G4PVPlacement *topPMTbase_phys					= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, topPMTbase_offsetZ),	topPMTbase_log, "topPMTbase_phys",	GXeVol_log, false, 0);
	G4PVPlacement *bottomPMTbase_phys				= new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, bottomPMTbase_offsetZ), bottomPMTbase_log, "bottomPMTbase_phys",	LXeVol_log, false, 0);

	// horizontal
	G4PVPlacement *NaI_Crystal_phys					= new G4PVPlacement(NaI_Crystal_T, NaI_Crystal_log, "NaI_Crystal_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaI_CrystalHousingTop_phys		= new G4PVPlacement(NaI_CrystalHousingTop_T, NaI_CrystalHousingTop_log, "NaI_CrystalHousingTop_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaI_CrystalHousingBot_phys		= new G4PVPlacement(NaI_CrystalHousingBot_T, NaI_CrystalHousingBot_log, "NaI_CrystalHousingBot_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaI_LightShieldBot_phys			= new G4PVPlacement(NaI_LightShieldBot_T, NaI_LightShieldBot_log, "NaI_LightShieldBot_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaI_LightShieldMid_phys			= new G4PVPlacement(NaI_LightShieldMid_T, NaI_LightShieldMid_log, "NaI_LightShieldMid_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaI_LightShieldTop_phys			= new G4PVPlacement(NaI_LightShieldTop_T, NaI_LightShieldTop_log, "NaI_LightShieldTop_phys", Laboratory_log, false, 0);
	//G4PVPlacement *NaI_LightShieldVac_phys		= new G4PVPlacement(NaI_LightShieldVac_T, NaI_LightShieldVac_log, "NaI_LightShieldVac_phys", Laboratory_log, false, 0);
	//G4PVPlacement *NaIcastle_phys					= new G4PVPlacement(0, G4ThreeVector(NaIcastle_offsetX, NaIcastle_offsetY, NaIcastle_offsetZ), NaIcastle_log, "NaIcastle_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIcastle_phys					= new G4PVPlacement(NaIcastle_T, NaIcastle_log, "NaIcastle_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannel_phys					= new G4PVPlacement(LeadChannel_T, LeadChannel_log, "LeadChannel_phys", Laboratory_log, false, 0);
	G4PVPlacement *Aperture_phys					= new G4PVPlacement(Aperture_T, Aperture_log, "Aperture_phys", Laboratory_log, false, 0);
	G4PVPlacement *ApertureFrame_phys				= new G4PVPlacement(ApertureFrame_T, ApertureFrame_log, "ApertureFrame_phys", Laboratory_log, false, 0);

	// NaI stand, with rotation
	G4PVPlacement *NaIstand_PlateTop_phys						= new G4PVPlacement(NaIstand_PlateTop_T, NaIstand_PlateTop_log, "NaIstand_PlateTop_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarOut_LeftFront_phys			= new G4PVPlacement(NaIstand_VertBar_LeftFront_T, NaIstand_VertBarOut_LeftFront_log, "NaIstand_VertBarOut_LeftFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarOut_LeftBack_phys			= new G4PVPlacement(NaIstand_VertBar_LeftBack_T, NaIstand_VertBarOut_LeftBack_log, "NaIstand_VertBarOut_LeftBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarOut_RightFront_phys			= new G4PVPlacement(NaIstand_VertBar_RightFront_T, NaIstand_VertBarOut_RightFront_log, "NaIstand_VertBarOut_RightFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarOut_RightBack_phys			= new G4PVPlacement(NaIstand_VertBar_RightBack_T, NaIstand_VertBarOut_RightBack_log, "NaIstand_VertBarOut_RightBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarIn_LeftFront_phys			= new G4PVPlacement(NaIstand_VertBar_LeftFront_T, NaIstand_VertBarIn_LeftFront_log, "NaIstand_VertBarIn_LeftFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarIn_LeftBack_phys				= new G4PVPlacement(NaIstand_VertBar_LeftBack_T, NaIstand_VertBarIn_LeftBack_log, "NaIstand_VertBarIn_LeftBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarIn_RightFront_phys			= new G4PVPlacement(NaIstand_VertBar_RightFront_T, NaIstand_VertBarIn_RightFront_log, "NaIstand_VertBarIn_RightFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_VertBarIn_RightBack_phys			= new G4PVPlacement(NaIstand_VertBar_RightBack_T, NaIstand_VertBarIn_RightBack_log, "NaIstand_VertBarIn_RightBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarOut_Left_phys				= new G4PVPlacement(NaIstand_HorzBar_Left_T, NaIstand_HorzBarOut_Left_log, "NaIstand_HorzBarOut_Left_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarOut_Right_phys				= new G4PVPlacement(NaIstand_HorzBar_Right_T, NaIstand_HorzBarOut_Right_log, "NaIstand_HorzBarOut_Right_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarOut_Front_phys				= new G4PVPlacement(NaIstand_HorzBar_Front_T, NaIstand_HorzBarOut_Front_log, "NaIstand_HorzBarOut_Front_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarOut_Back_phys				= new G4PVPlacement(NaIstand_HorzBar_Back_T, NaIstand_HorzBarOut_Back_log, "NaIstand_HorzBarOut_Back_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarIn_Left_phys					= new G4PVPlacement(NaIstand_HorzBar_Left_T, NaIstand_HorzBarIn_Left_log, "NaIstand_HorzBarIn_Left_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarIn_Right_phys				= new G4PVPlacement(NaIstand_HorzBar_Right_T, NaIstand_HorzBarIn_Right_log, "NaIstand_HorzBarIn_Right_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarIn_Front_phys				= new G4PVPlacement(NaIstand_HorzBar_Front_T, NaIstand_HorzBarIn_Front_log, "NaIstand_HorzBarIn_Front_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_HorzBarIn_Back_phys					= new G4PVPlacement(NaIstand_HorzBar_Back_T, NaIstand_HorzBarIn_Back_log, "NaIstand_HorzBarIn_Back_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_PlateBot_phys						= new G4PVPlacement(NaIstand_PlateBot_T, NaIstand_PlateBot_log, "NaIstand_PlateBot_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_CradleBarLeft_phys					= new G4PVPlacement(NaIstand_CradleBarLeft_T, NaIstand_CradleBarLeft_log, "NaIstand_CradleBarLeft_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_CradleBarRight_phys					= new G4PVPlacement(NaIstand_CradleBarRight_T, NaIstand_CradleBarRight_log, "NaIstand_CradleBarRight_phys", Laboratory_log, false, 0);
	G4PVPlacement *NaIstand_CradlePanel_phys					= new G4PVPlacement(NaIstand_CradlePanel_T, NaIstand_CradlePanel_log, "NaIstand_CradlePanel_phys", Laboratory_log, false, 0);

	// Lead Channel stand, with rotation; replicated from NaIstand
	G4PVPlacement *LeadChannelStand_PlateTop_phys				= new G4PVPlacement(LeadChannelStand_PlateTop_T, LeadChannelStand_PlateTop_log, "LeadChannelStand_PlateTop_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarOut_LeftFront_phys	= new G4PVPlacement(LeadChannelStand_VertBar_LeftFront_T, LeadChannelStand_VertBarOut_LeftFront_log, "LeadChannelStand_VertBarOut_LeftFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarOut_LeftBack_phys	= new G4PVPlacement(LeadChannelStand_VertBar_LeftBack_T, LeadChannelStand_VertBarOut_LeftBack_log, "LeadChannelStand_VertBarOut_LeftBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarOut_RightFront_phys	= new G4PVPlacement(LeadChannelStand_VertBar_RightFront_T, LeadChannelStand_VertBarOut_RightFront_log, "LeadChannelStand_VertBarOut_RightFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarOut_RightBack_phys	= new G4PVPlacement(LeadChannelStand_VertBar_RightBack_T, LeadChannelStand_VertBarOut_RightBack_log, "LeadChannelStand_VertBarOut_RightBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarIn_LeftFront_phys	= new G4PVPlacement(LeadChannelStand_VertBar_LeftFront_T, LeadChannelStand_VertBarIn_LeftFront_log, "LeadChannelStand_VertBarIn_LeftFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarIn_LeftBack_phys		= new G4PVPlacement(LeadChannelStand_VertBar_LeftBack_T, LeadChannelStand_VertBarIn_LeftBack_log, "LeadChannelStand_VertBarIn_LeftBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarIn_RightFront_phys	= new G4PVPlacement(LeadChannelStand_VertBar_RightFront_T, LeadChannelStand_VertBarIn_RightFront_log, "LeadChannelStand_VertBarIn_RightFront_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_VertBarIn_RightBack_phys	= new G4PVPlacement(LeadChannelStand_VertBar_RightBack_T, LeadChannelStand_VertBarIn_RightBack_log, "LeadChannelStand_VertBarIn_RightBack_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarOut_Left_phys		= new G4PVPlacement(LeadChannelStand_HorzBar_Left_T, LeadChannelStand_HorzBarOut_Left_log, "LeadChannelStand_HorzBarOut_Left_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarOut_Right_phys		= new G4PVPlacement(LeadChannelStand_HorzBar_Right_T, LeadChannelStand_HorzBarOut_Right_log, "LeadChannelStand_HorzBarOut_Right_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarOut_Front_phys		= new G4PVPlacement(LeadChannelStand_HorzBar_Front_T, LeadChannelStand_HorzBarOut_Front_log, "LeadChannelStand_HorzBarOut_Front_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarOut_Back_phys		= new G4PVPlacement(LeadChannelStand_HorzBar_Back_T, LeadChannelStand_HorzBarOut_Back_log, "LeadChannelStand_HorzBarOut_Back_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarIn_Left_phys			= new G4PVPlacement(LeadChannelStand_HorzBar_Left_T, LeadChannelStand_HorzBarIn_Left_log, "LeadChannelStand_HorzBarIn_Left_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarIn_Right_phys		= new G4PVPlacement(LeadChannelStand_HorzBar_Right_T, LeadChannelStand_HorzBarIn_Right_log, "LeadChannelStand_HorzBarIn_Right_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarIn_Front_phys		= new G4PVPlacement(LeadChannelStand_HorzBar_Front_T, LeadChannelStand_HorzBarIn_Front_log, "LeadChannelStand_HorzBarIn_Front_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_HorzBarIn_Back_phys			= new G4PVPlacement(LeadChannelStand_HorzBar_Back_T, LeadChannelStand_HorzBarIn_Back_log, "LeadChannelStand_HorzBarIn_Back_phys", Laboratory_log, false, 0);
	G4PVPlacement *LeadChannelStand_PlateBot_phys				= new G4PVPlacement(LeadChannelStand_PlateBot_T, LeadChannelStand_PlateBot_log, "LeadChannelStand_PlateBot_phys", Laboratory_log, false, 0);

	// collimator
	G4PVPlacement *Collimator_phys						= new G4PVPlacement(Collimator_T, Collimator_log, "Collimator_phys", Laboratory_log, false, 0);

	////////////////////////////////////////////////////////////////////////
	// VISUAL ATTRIBUTES
	////////////////////////////////////////////////////////////////////////
	Laboratory_log			->SetVisAttributes(G4VisAttributes::Invisible);
	CryostatVacuum_log		->SetVisAttributes(G4VisAttributes::Invisible);
	topPMTinterior_log		->SetVisAttributes(G4VisAttributes::Invisible);
	bottomPMTinterior_log	->SetVisAttributes(G4VisAttributes::Invisible);
	//NaI_LightShieldVac_log	->SetVisAttributes(G4VisAttributes::Invisible);

	G4VisAttributes *Al_vis = new G4VisAttributes(grey);
		XuStand_HorzBar_TopLeft_log					->SetVisAttributes(Al_vis);
		XuStand_HorzBar_TopRight_log				->SetVisAttributes(Al_vis);
		XuStand_HorzBar_TopFront_log				->SetVisAttributes(Al_vis);
		XuStand_HorzBar_TopBack_log					->SetVisAttributes(Al_vis);
		XuStand_WingFront_log						->SetVisAttributes(Al_vis);
		XuStand_WingBack_log						->SetVisAttributes(Al_vis);
		XuStand_HorzBar_MidLeft_log					->SetVisAttributes(Al_vis);
		XuStand_HorzBar_MidRight_log				->SetVisAttributes(Al_vis);
		XuStand_HorzBar_MidFront_log				->SetVisAttributes(Al_vis);
		XuStand_HorzBar_MidBack_log					->SetVisAttributes(Al_vis);
		XuStand_VertBar_LeftFront_log				->SetVisAttributes(Al_vis);
		XuStand_VertBar_LeftBack_log				->SetVisAttributes(Al_vis);
		XuStand_VertBar_RightFront_log				->SetVisAttributes(Al_vis);
		XuStand_VertBar_RightBack_log				->SetVisAttributes(Al_vis);
		XuStand_HorzBar_BotLeft_log					->SetVisAttributes(Al_vis);
		XuStand_HorzBar_BotRight_log				->SetVisAttributes(Al_vis);
		XuStand_HorzBar_BotFront_log				->SetVisAttributes(Al_vis);
		XuStand_HorzBar_BotBack_log					->SetVisAttributes(Al_vis);
		XuStand_BottomPlate_log						->SetVisAttributes(Al_vis);
		XuStand_CylinderLeftFront_log				->SetVisAttributes(Al_vis);
		XuStand_CylinderLeftBack_log				->SetVisAttributes(Al_vis);
		XuStand_CylinderRightFront_log				->SetVisAttributes(Al_vis);
		XuStand_CylinderRightBack_log				->SetVisAttributes(Al_vis);
		XuStand_CradleBarLeftFront_log				->SetVisAttributes(Al_vis);
		XuStand_CradleBarLeftBack_log				->SetVisAttributes(Al_vis);
		XuStand_CradleBarRightFront_log				->SetVisAttributes(Al_vis);
		XuStand_CradleBarRightBack_log				->SetVisAttributes(Al_vis);
		XuStand_CradlePanelLong_log					->SetVisAttributes(Al_vis);
		XuStand_CradlePanelShort_log				->SetVisAttributes(Al_vis);
		XuStand_Goniometer_log						->SetVisAttributes(Al_vis);
		NaIstand_PlateTop_log						->SetVisAttributes(Al_vis);
		NaIstand_VertBarOut_LeftFront_log			->SetVisAttributes(Al_vis);
		NaIstand_VertBarOut_LeftBack_log			->SetVisAttributes(Al_vis);
		NaIstand_VertBarOut_RightFront_log			->SetVisAttributes(Al_vis);
		NaIstand_VertBarOut_RightBack_log			->SetVisAttributes(Al_vis);
		NaIstand_VertBarIn_LeftFront_log			->SetVisAttributes(Al_vis);
		NaIstand_VertBarIn_LeftBack_log				->SetVisAttributes(Al_vis);
		NaIstand_VertBarIn_RightFront_log			->SetVisAttributes(Al_vis);
		NaIstand_VertBarIn_RightBack_log			->SetVisAttributes(Al_vis);
		NaIstand_HorzBarOut_Left_log				->SetVisAttributes(Al_vis);
		NaIstand_HorzBarOut_Right_log				->SetVisAttributes(Al_vis);		
		NaIstand_HorzBarOut_Front_log				->SetVisAttributes(Al_vis);	
		NaIstand_HorzBarOut_Back_log				->SetVisAttributes(Al_vis);
		NaIstand_HorzBarIn_Left_log					->SetVisAttributes(Al_vis);
		NaIstand_HorzBarIn_Right_log				->SetVisAttributes(Al_vis);
		NaIstand_HorzBarIn_Front_log				->SetVisAttributes(Al_vis);
		NaIstand_HorzBarIn_Back_log					->SetVisAttributes(Al_vis);
		NaIstand_PlateBot_log						->SetVisAttributes(Al_vis);
		NaIstand_CradleBarLeft_log					->SetVisAttributes(Al_vis);
		NaIstand_CradleBarRight_log					->SetVisAttributes(Al_vis);
		NaIstand_CradlePanel_log					->SetVisAttributes(Al_vis);
		LeadChannelStand_PlateTop_log				->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarOut_LeftFront_log	->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarOut_LeftBack_log	->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarOut_RightFront_log	->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarOut_RightBack_log	->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarIn_LeftFront_log	->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarIn_LeftBack_log		->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarIn_RightFront_log	->SetVisAttributes(Al_vis);
		LeadChannelStand_VertBarIn_RightBack_log	->SetVisAttributes(Al_vis);
		LeadChannelStand_HorzBarOut_Left_log		->SetVisAttributes(Al_vis);
		LeadChannelStand_HorzBarOut_Right_log		->SetVisAttributes(Al_vis);		
		LeadChannelStand_HorzBarOut_Front_log		->SetVisAttributes(Al_vis);	
		LeadChannelStand_HorzBarOut_Back_log		->SetVisAttributes(Al_vis);
		LeadChannelStand_HorzBarIn_Left_log			->SetVisAttributes(Al_vis);
		LeadChannelStand_HorzBarIn_Right_log		->SetVisAttributes(Al_vis);
		LeadChannelStand_HorzBarIn_Front_log		->SetVisAttributes(Al_vis);
		LeadChannelStand_HorzBarIn_Back_log			->SetVisAttributes(Al_vis);
		LeadChannelStand_PlateBot_log				->SetVisAttributes(Al_vis);
		AlCanTube_log								->SetVisAttributes(Al_vis);
		AlCanBottom_log								->SetVisAttributes(Al_vis);
		ApertureFrame_log							->SetVisAttributes(Al_vis);

	G4VisAttributes *LXe_vis = new G4VisAttributes(xlred);
	LXe_vis	->SetVisibility(true);
	LXe_vis	->SetForceSolid(false);
		LXeVol_log		->SetVisAttributes(LXe_vis);
		GXeVol_log		->SetVisAttributes(LXe_vis);
	
	G4VisAttributes	*Teflon_vis	= new G4VisAttributes(magenta);
	Teflon_vis	->SetVisibility(true);
	Teflon_vis	->SetForceSolid(false);
		Washer1_log						->SetVisAttributes(Teflon_vis);
		Washer2_log						->SetVisAttributes(Teflon_vis);
		Washer3_log						->SetVisAttributes(Teflon_vis);
		Washer4_log						->SetVisAttributes(Teflon_vis);
		Support1_log					->SetVisAttributes(Teflon_vis);
		Support2_log					->SetVisAttributes(Teflon_vis);
		Support3_log					->SetVisAttributes(Teflon_vis);
		Support4_log					->SetVisAttributes(Teflon_vis);
		TopClamp_log					->SetVisAttributes(Teflon_vis);
		TopHolderUpperTube_log			->SetVisAttributes(Teflon_vis);
		TopHolderLowerTube_log			->SetVisAttributes(Teflon_vis);
		TopHolderSquare_log				->SetVisAttributes(Teflon_vis);
		ExtractionSpacerUpperHalf_log	->SetVisAttributes(Teflon_vis);
		ExtractionSpacerLowerHalf_log	->SetVisAttributes(Teflon_vis);
		DriftSpacerUpperSquare_log		->SetVisAttributes(Teflon_vis);
		DriftSpacerUpperTube_log		->SetVisAttributes(Teflon_vis);
		DriftSpacerMiddleTube_log		->SetVisAttributes(Teflon_vis);
		DriftSpacerLowerTube_log		->SetVisAttributes(Teflon_vis);
		DriftSpacerLowerSquare_log		->SetVisAttributes(Teflon_vis);
		BottomHolderUpperTube_log		->SetVisAttributes(Teflon_vis);
		BottomHolderMiddleTube_log		->SetVisAttributes(Teflon_vis);
		BottomHolderLowerTube_log		->SetVisAttributes(Teflon_vis);
		Filler_log						->SetVisAttributes(Teflon_vis);
		BottomClampUpperPart_log		->SetVisAttributes(Teflon_vis);
		BottomClampMiddlePart_log		->SetVisAttributes(Teflon_vis);
		BottomClampLowerPart_log		->SetVisAttributes(Teflon_vis);
		topPMTbase_log					->SetVisAttributes(Teflon_vis);
		bottomPMTbase_log				->SetVisAttributes(Teflon_vis);
		
	G4VisAttributes *Cryostat_vis = new G4VisAttributes(xlblue);
	Cryostat_vis	->SetVisibility(true);
	Cryostat_vis	->SetForceSolid(false);
		OuterCanTube_log		->SetVisAttributes(Cryostat_vis);
		OuterCanBottom_log		->SetVisAttributes(Cryostat_vis);
		OuterCanLowerFlange_log	->SetVisAttributes(Cryostat_vis);
		OuterCanUpperFlange_log	->SetVisAttributes(Cryostat_vis);
		InnerCanTube_log		->SetVisAttributes(Cryostat_vis);
		InnerCanBottom_log		->SetVisAttributes(Cryostat_vis);
		InnerCanLowerFlange_log	->SetVisAttributes(Cryostat_vis);
		InnerCanUpperFlange_log	->SetVisAttributes(Cryostat_vis);
		
		Pipe1in_log		->SetVisAttributes(Cryostat_vis);
		Pipe2in_log		->SetVisAttributes(Cryostat_vis);
		Pipe3in_log		->SetVisAttributes(Cryostat_vis);
		Pipe4in_log		->SetVisAttributes(Cryostat_vis);
		Pipe5in_log		->SetVisAttributes(Cryostat_vis);
		Pipe6in_log		->SetVisAttributes(Cryostat_vis);
		Pipe7in_log		->SetVisAttributes(Cryostat_vis);
		Pipe8in_log		->SetVisAttributes(Cryostat_vis);
		Pipe9in_log		->SetVisAttributes(Cryostat_vis);
		Pipe10in_log	->SetVisAttributes(Cryostat_vis);
		Pipe11in_log	->SetVisAttributes(Cryostat_vis);

		Pipe1out_log	->SetVisAttributes(Cryostat_vis);
		Pipe2out_log	->SetVisAttributes(Cryostat_vis);
		Pipe3out_log	->SetVisAttributes(Cryostat_vis);
		Pipe4out_log	->SetVisAttributes(Cryostat_vis);
		Pipe5out_log	->SetVisAttributes(Cryostat_vis);
		Pipe6out_log	->SetVisAttributes(Cryostat_vis);
		Pipe7out_log	->SetVisAttributes(Cryostat_vis);
		Pipe8out_log	->SetVisAttributes(Cryostat_vis);
		Pipe9out_log	->SetVisAttributes(Cryostat_vis);
		Pipe10out_log	->SetVisAttributes(Cryostat_vis);
		Pipe11out_log	->SetVisAttributes(Cryostat_vis);

		Pipe1flange_log	->SetVisAttributes(Cryostat_vis);
		Pipe2flange_log	->SetVisAttributes(Cryostat_vis);
		Pipe3flange_log	->SetVisAttributes(Cryostat_vis);
		Pipe4flange_log	->SetVisAttributes(Cryostat_vis);
		Pipe5flange1_log->SetVisAttributes(Cryostat_vis);
		Pipe5flange2_log->SetVisAttributes(Cryostat_vis);
		Pipe6flange_log	->SetVisAttributes(Cryostat_vis);
		Pipe8flange_log	->SetVisAttributes(Cryostat_vis);
		Pipe9flange_log	->SetVisAttributes(Cryostat_vis);
		Pipe10flange_log->SetVisAttributes(Cryostat_vis);
		Pipe11flange_log->SetVisAttributes(Cryostat_vis);

		DewarOutTube_log->SetVisAttributes(Cryostat_vis);
		DewarOutBot_log	->SetVisAttributes(Cryostat_vis);
		DewarInTube_log	->SetVisAttributes(Cryostat_vis);
		DewarInBot_log	->SetVisAttributes(Cryostat_vis);
		DewarTop_log	->SetVisAttributes(Cryostat_vis);

	G4VisAttributes	*SteelHolder_vis = new G4VisAttributes(blue);
	SteelHolder_vis	->SetVisibility(true);
	SteelHolder_vis	->SetForceSolid(false);	
		SteelHolderTube_log			->SetVisAttributes(SteelHolder_vis);
		SteelHolderUpperSquare_log	->SetVisAttributes(SteelHolder_vis);
		SteelHolderLowerSquare_log	->SetVisAttributes(SteelHolder_vis);
		
	G4VisAttributes	*PMTcasing_vis = new G4VisAttributes(cyan);
	PMTcasing_vis	->SetVisibility(true);
	PMTcasing_vis	->SetForceSolid(false);
		topPMTcasing_log	->SetVisAttributes(PMTcasing_vis);	
		bottomPMTcasing_log	->SetVisAttributes(PMTcasing_vis);	

	G4VisAttributes	*PMTinterior_vis = new G4VisAttributes(grey);
	PMTinterior_vis	->SetVisibility(true);
	PMTinterior_vis	->SetForceSolid(false);
		//topPMTinterior_log		->SetVisAttributes(PMTinterior_vis);	
		//bottomPMTinterior_log	->SetVisAttributes(PMTinterior_vis);	

	G4VisAttributes	* PMTwindow_vis	= new G4VisAttributes(grey);
	PMTwindow_vis	->SetVisibility(true);
	PMTwindow_vis	->SetForceSolid(false);
		topPMTwindow_log	->SetVisAttributes(PMTwindow_vis);
		bottomPMTwindow_log	->SetVisAttributes(PMTwindow_vis);

	G4VisAttributes	*PMTcathode_vis	= new G4VisAttributes(xlgreen);
	PMTcathode_vis	->SetVisibility(true);
	PMTcathode_vis	->SetForceSolid(false);
		topPMTcathode_log		->SetVisAttributes(PMTcathode_vis);
		bottomPMTcathode_log	->SetVisAttributes(PMTcathode_vis);
		
	G4VisAttributes	*Grids_vis = new G4VisAttributes(blue);
	Grids_vis	->SetVisibility(true);
	Grids_vis	->SetForceSolid(false);
		AnodeGrid_log	->SetVisAttributes(Grids_vis);
		GateGrid_log	->SetVisAttributes(Grids_vis);
		CathodeGrid_log	->SetVisAttributes(Grids_vis);

	G4VisAttributes	*Cu_vis = new G4VisAttributes(yellow);
	Cu_vis	->SetVisibility(true);
	Cu_vis	->SetForceSolid(false);
		CopperFinger_log	->SetVisAttributes(Cu_vis);

	G4VisAttributes	*NaI_Crystal_vis = new G4VisAttributes(xlred);
	NaI_Crystal_vis	->SetVisibility(true);
	NaI_Crystal_vis	->SetForceSolid(false);
		NaI_Crystal_log	->SetVisAttributes(NaI_Crystal_vis);

	G4VisAttributes	*NaI_CrystalHousing_vis = new G4VisAttributes(lblue);
	NaI_CrystalHousing_vis	->SetVisibility(true);
	NaI_CrystalHousing_vis	->SetForceSolid(false);
		NaI_CrystalHousingTop_log	->SetVisAttributes(NaI_CrystalHousing_vis);
		NaI_CrystalHousingBot_log	->SetVisAttributes(NaI_CrystalHousing_vis);

	G4VisAttributes	*NaI_LightShield_vis = new G4VisAttributes(blue);
	NaI_LightShield_vis	->SetVisibility(true);
	NaI_LightShield_vis	->SetForceSolid(false);
		NaI_LightShieldTop_log	->SetVisAttributes(NaI_LightShield_vis);
		NaI_LightShieldBot_log	->SetVisAttributes(NaI_LightShield_vis);
		NaI_LightShieldMid_log	->SetVisAttributes(NaI_LightShield_vis);

	G4VisAttributes	*NaI_PMTcasing_vis = new G4VisAttributes(lgreen);
	NaI_PMTcasing_vis	->SetVisibility(true);
	NaI_PMTcasing_vis	->SetForceSolid(false);
		//NaI_PMTcasing_log	->SetVisAttributes(NaI_PMTcasing_vis);

	G4VisAttributes	*Collimator_vis = new G4VisAttributes(lgreen);
		Collimator_vis	->SetVisibility(true);
		Collimator_vis	->SetForceSolid(false);
		Collimator_log	->SetVisAttributes(Collimator_vis);

	G4VisAttributes	*NaIcastle_vis = new G4VisAttributes(lgreen);
		NaIcastle_vis	->SetVisibility(true);
		NaIcastle_vis	->SetForceSolid(false);
		NaIcastle_log	->SetVisAttributes(NaIcastle_vis);

	G4VisAttributes	*LeadChannel_vis = new G4VisAttributes(xlgreen);
		LeadChannel_vis	->SetVisibility(true);
		LeadChannel_vis	->SetForceSolid(false);
		LeadChannel_log	->SetVisAttributes(LeadChannel_vis);
		Aperture_log	->SetVisAttributes(LeadChannel_vis);

		
//--------------	DETECTOR DESCRIPTION	------------
	G4double Mass_Laboratory				= Laboratory_log				->GetMass(false, false, 0)/kg;	
	G4double Mass_OuterCanTube				= OuterCanTube_log				->GetMass(false, false, 0)/kg;		
	G4double Mass_OuterCanBottom			= OuterCanBottom_log			->GetMass(false, false, 0)/kg;		
	G4double Mass_OuterCanLowerFlange		= OuterCanLowerFlange_log		->GetMass(false, false, 0)/kg;		
	G4double Mass_OuterCanUpperFlange		= OuterCanUpperFlange_log		->GetMass(false, false, 0)/kg;		
	G4double Mass_InnerCanTube				= InnerCanTube_log				->GetMass(false, false, 0)/kg;
	G4double Mass_InnerCanBottom			= InnerCanBottom_log			->GetMass(false, false, 0)/kg;
	G4double Mass_InnerCanLowerFlange		= InnerCanLowerFlange_log		->GetMass(false, false, 0)/kg;
	G4double Mass_InnerCanUpperFlange		= InnerCanUpperFlange_log		->GetMass(false, false, 0)/kg;
	G4double Mass_SteelHolderUpperSquare	= SteelHolderUpperSquare_log	->GetMass(false, false, 0)/kg;		
	G4double Mass_SteelHolderLowerSquare	= SteelHolderLowerSquare_log	->GetMass(false, false, 0)/kg;		
	G4double Mass_SteelHolderTube			= SteelHolderTube_log			->GetMass(false, false, 0)/kg;		
	G4double Mass_Support					= Support1_log					->GetMass(false, false, 0)/kg;
	G4double Mass_Washer					= Washer1_log					->GetMass(false, false, 0)/kg;
	G4double Mass_TopClamp					= TopClamp_log					->GetMass(false, false, 0)/kg;
	G4double Mass_TopHolderUpperTube		= TopHolderUpperTube_log		->GetMass(false, false, 0)/kg;
	G4double Mass_TopHolderLowerTube		= TopHolderLowerTube_log		->GetMass(false, false, 0)/kg;
	G4double Mass_TopHolderSquare			= TopHolderSquare_log			->GetMass(false, false, 0)/kg;
	G4double Mass_AnodeGrid					= AnodeGrid_log					->GetMass(false, false, 0)/kg;
	G4double Mass_ExtractionSpacerUpperHalf	= ExtractionSpacerUpperHalf_log	->GetMass(false, false, 0)/kg;
	G4double Mass_ExtractionSpacerLowerHalf	= ExtractionSpacerLowerHalf_log	->GetMass(false, false, 0)/kg;
	G4double Mass_GateGrid					= GateGrid_log					->GetMass(false, false, 0)/kg;	
	G4double Mass_DriftSpacerUpperSquare	= DriftSpacerUpperSquare_log	->GetMass(false, false, 0)/kg;
	G4double Mass_DriftSpacerUpperTube		= DriftSpacerUpperTube_log		->GetMass(false, false, 0)/kg;
	G4double Mass_DriftSpacerMiddleTube		= DriftSpacerMiddleTube_log		->GetMass(false, false, 0)/kg;
	G4double Mass_DriftSpacerLowerTube		= DriftSpacerLowerTube_log		->GetMass(false, false, 0)/kg;
	G4double Mass_DriftSpacerLowerSquare	= DriftSpacerLowerSquare_log	->GetMass(false, false, 0)/kg;
	G4double Mass_BottomHolderUpperTube		= BottomHolderUpperTube_log		->GetMass(false, false, 0)/kg;
	G4double Mass_BottomHolderMiddleTube	= BottomHolderMiddleTube_log	->GetMass(false, false, 0)/kg;
	G4double Mass_BottomHolderLowerTube		= BottomHolderLowerTube_log		->GetMass(false, false, 0)/kg;
	G4double Mass_BottomClampUpperPart		= BottomClampUpperPart_log		->GetMass(false, false, 0)/kg;
	G4double Mass_BottomClampMiddlePart		= BottomClampMiddlePart_log		->GetMass(false, false, 0)/kg;
	G4double Mass_BottomClampLowerPart		= BottomClampLowerPart_log		->GetMass(false, false, 0)/kg;
	G4double Mass_Filler					= Filler_log					->GetMass(false, false, 0)/kg;
	G4double Mass_topPMTbase				= topPMTbase_log				->GetMass(false, false, 0)/kg;
	G4double Mass_topPMTcasing				= topPMTcasing_log				->GetMass(false, false, 0)/kg;
	G4double Mass_topPMTinterior			= topPMTinterior_log			->GetMass(false, false, 0)/kg;
	G4double Mass_topPMTwindow				= topPMTwindow_log				->GetMass(false, false, 0)/kg;
	G4double Mass_topPMTcathode				= topPMTcathode_log				->GetMass(false, false, 0)/kg;
	G4double Mass_CathodeGrid				= CathodeGrid_log				->GetMass(false, false, 0)/kg;
	G4double Mass_GXe						= GXeVol_log					->GetMass(false, false, 0)/kg;
	G4double Mass_LXe						= LXeVol_log					->GetMass(false, false, 0)/kg;
	G4double Mass_Target					= LXeTarget_log					->GetMass(false, false, 0)/kg;
	G4double Mass_Gate						= LXeGate_log					->GetMass(false, false, 0)/kg;
	//G4double Mass_Cathode					= LXeCathode_log				->GetMass(false, false, 0)/kg;
	G4double Mass_TargetTotal				= Mass_Target + Mass_Gate;

	// feedthroughs
	G4double Mass_Pipes = 0;
			 Mass_Pipes += Pipe1in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe2in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe3in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe4in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe5in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe6in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe7in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe8in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe9in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe10in_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe11in_log		->GetMass(false, false, 0)/kg;

			 Mass_Pipes += Pipe1out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe2out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe3out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe4out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe5out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe6out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe7out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe8out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe9out_log		->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe10out_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe11out_log	->GetMass(false, false, 0)/kg;

			 Mass_Pipes += Pipe1flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe2flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe3flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe4flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe5flange1_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe5flange2_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe6flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe8flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe9flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe10flange_log	->GetMass(false, false, 0)/kg;
			 Mass_Pipes += Pipe11flange_log	->GetMass(false, false, 0)/kg;

	// Xuerich stand
	G4double Mass_XuStand = 0;
			 Mass_XuStand += XuStand_HorzBar_TopLeft_log	->GetMass(false, false, 0)/kg;
			 Mass_XuStand += XuStand_HorzBar_TopRight_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_TopFront_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_TopBack_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_WingFront_log			->GetMass(false, false, 0)/kg;			
			 Mass_XuStand += XuStand_WingBack_log			->GetMass(false, false, 0)/kg;			
			 Mass_XuStand += XuStand_HorzBar_MidLeft_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_MidRight_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_MidFront_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_MidBack_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_VertBar_LeftFront_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_VertBar_LeftBack_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_VertBar_RightFront_log	->GetMass(false, false, 0)/kg;
			 Mass_XuStand += XuStand_VertBar_RightBack_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_BotLeft_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_BotRight_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_BotFront_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_HorzBar_BotBack_log	->GetMass(false, false, 0)/kg;	
			 Mass_XuStand += XuStand_BottomPlate_log		->GetMass(false, false, 0)/kg;		
	// copper finger
	G4double Mass_CopperFinger = CopperFinger_log	->GetMass(false, false, 0)/kg;
	// golometer extension
	G4double Mass_XuStandExtension = 0;
			 Mass_XuStandExtension += XuStand_CylinderLeftFront_log		->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CylinderLeftBack_log		->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CylinderRightFront_log	->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CylinderRightBack_log		->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CradleBarLeftFront_log	->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CradleBarLeftBack_log		->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CradleBarRightFront_log	->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CradleBarRightBack_log	->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CradlePanelLong_log		->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_CradlePanelShort_log		->GetMass(false, false, 0)/kg;
			 Mass_XuStandExtension += XuStand_Goniometer_log			->GetMass(false, false, 0)/kg;
	G4double Mass_NaIstand = 0;		
			 Mass_NaIstand += NaIstand_PlateTop_log				->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarOut_LeftFront_log	->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarOut_LeftBack_log	->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarOut_RightFront_log->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarOut_RightBack_log	->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarIn_LeftFront_log	->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarIn_LeftBack_log	->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarIn_RightFront_log	->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_VertBarIn_RightBack_log	->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarOut_Left_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarOut_Right_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarOut_Front_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarOut_Back_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarIn_Left_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarIn_Right_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarIn_Front_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_HorzBarIn_Back_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_PlateBot_log				->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_CradleBarLeft_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_CradleBarRight_log		->GetMass(false, false, 0)/kg;
			 Mass_NaIstand += NaIstand_CradlePanel_log			->GetMass(false, false, 0)/kg;
	
	G4cout	<< "|   DETECTOR PARTS (AS CODED) AND WEIGHT: "													<< G4endl
			<< "|   Laboratory                       "	<<	Mass_Laboratory					<< " kg"		<< G4endl
			<< "|---Al STANDS------------------------"														<< G4endl
			<< "|   Standard Xuerich stand           "	<<	Mass_XuStand					<< " kg"		<< G4endl
			<< "|   Goniometer extension for Xuerich "	<<	Mass_XuStandExtension			<< " kg"		<< G4endl
			<< "|   Stand for NaI detector           "	<<	Mass_NaIstand					<< " kg"		<< G4endl
			<< "|---Copper Finger--------------------"														<< G4endl
			<< "|   Copper Finger                    "	<<	Mass_CopperFinger				<< " kg"		<< G4endl
			<< "|---STEEL----------------------------"														<< G4endl
			<< "|   Outer Can Tube                   "	<<	Mass_OuterCanTube				<< " kg"		<< G4endl
			<< "|   Outer Can Bottom                 "	<<	Mass_OuterCanBottom				<< " kg"		<< G4endl
			<< "|   Outer Can Upper Flange           "	<<	Mass_OuterCanLowerFlange		<< " kg"		<< G4endl
			<< "|   Outer Can Lower Flange           "	<<	Mass_OuterCanUpperFlange		<< " kg"		<< G4endl
			<< "|   Inner Can Tube                   "	<<	Mass_InnerCanTube				<< " kg"		<< G4endl
			<< "|   Inner Can Bottom                 "	<<	Mass_InnerCanBottom				<< " kg"		<< G4endl
			<< "|   Inner Can LowerFlange            "	<<	Mass_InnerCanLowerFlange		<< " kg"		<< G4endl
			<< "|   Inner Can Upper Flange           "	<<	Mass_InnerCanUpperFlange		<< " kg"		<< G4endl
			<< "|   Steel Holder Tube                "	<<	Mass_SteelHolderTube			<< " kg"		<< G4endl
			<< "|   Steel Holder Upper Square        "	<<	Mass_SteelHolderUpperSquare		<< " kg"		<< G4endl
			<< "|   Steel Holder Lower Square        "	<<	Mass_SteelHolderLowerSquare		<< " kg"		<< G4endl
			<< "|   Feedthroughs                     "	<<	Mass_Pipes						<< " kg"		<< G4endl
			<< "|---TEFLON---------------------------"														<< G4endl
			<< "|   Support (x 4)                    "	<<	Mass_Support					<< " kg"		<< G4endl
			<< "|   Washer (x 4)                     "	<<  Mass_Washer						<< " kg"		<< G4endl
			<< "|   Top PMT Clamp                    "	<<	Mass_TopClamp					<< " kg"		<< G4endl
			<< "|   Top PMT Holder Upper Tube        "	<<	Mass_TopHolderUpperTube			<< " kg"		<< G4endl
			<< "|   Top PMT Holder Lower Tube        "	<<	Mass_TopHolderLowerTube			<< " kg"		<< G4endl
			<< "|   Top PMT Holder Square            "	<<	Mass_TopHolderSquare			<< " kg"		<< G4endl
			<< "|   Anode Grid Frame                 "	<<	Mass_AnodeGrid					<< " kg"		<< G4endl
			<< "|   Extraction Spacer Upper Half     "	<<	Mass_ExtractionSpacerUpperHalf	<< " kg"		<< G4endl
			<< "|   Extraction Spacer Lower Half     "	<<	Mass_ExtractionSpacerLowerHalf	<< " kg"		<< G4endl
			<< "|   Gate Grid Frame                  "	<<	Mass_GateGrid					<< " kg"		<< G4endl
			<< "|   Drift Spacer Upper Square        "	<<	Mass_DriftSpacerUpperSquare		<< " kg"		<< G4endl
			<< "|   Drift Spacer Upper Tube          "	<<	Mass_DriftSpacerUpperTube		<< " kg"		<< G4endl
			<< "|   Drift Spacer Middle Tube         "	<<	Mass_DriftSpacerMiddleTube		<< " kg"		<< G4endl
			<< "|   Drift Spacer Lower Tube          "	<<	Mass_DriftSpacerLowerTube		<< " kg"		<< G4endl
			<< "|   Drift Spacer Lower Square        "	<<	Mass_DriftSpacerLowerSquare		<< " kg"		<< G4endl
			<< "|   Bottom PMT Holder Upper Tube     "	<<	Mass_BottomHolderUpperTube		<< " kg"		<< G4endl
			<< "|   Bottom PMT Holder Middle Tube    "	<<	Mass_BottomHolderMiddleTube		<< " kg"		<< G4endl
			<< "|   Bottom PMT Holder Lower Tube     "	<<	Mass_BottomHolderLowerTube		<< " kg"		<< G4endl
			<< "|   Bottom PMT Clamp Upper Part      "	<<	Mass_BottomClampUpperPart		<< " kg"		<< G4endl
			<< "|   Bottom PMT Clamp Middle Part     "	<<	Mass_BottomClampMiddlePart		<< " kg"		<< G4endl
			<< "|   Bottom PMT Clamp Lower Part      "	<<	Mass_BottomClampLowerPart		<< " kg"		<< G4endl
			<< "|   Filler                           "	<<	Mass_Filler						<< " kg"		<< G4endl
			<< "|   Grid Frame                       "	<<	Mass_CathodeGrid				<< " kg"		<< G4endl
			<< "|---PMT------------------------------"														<< G4endl
			<< "|   PMT base                         "	<<	Mass_topPMTbase					<< " kg"		<< G4endl
			<< "|   PMT casing                       "	<<	Mass_topPMTcasing				<< " kg"		<< G4endl
			<< "|   PMT interior                     "	<<	Mass_topPMTinterior				<< " kg"		<< G4endl
			<< "|   PMT window                       "	<<	Mass_topPMTwindow				<< " kg"		<< G4endl
			<< "|   PMT cathode                      "	<<	Mass_topPMTcathode				<< " kg"		<< G4endl
			<< "|---XENON----------------------------"														<< G4endl
			<< "|   GXe                              "	<<	Mass_GXe						<< " kg"		<< G4endl
			<< "|   LXe                              "	<<	Mass_LXe						<< " kg"		<< G4endl
			<< "|   LXe_Target                       "	<<	Mass_Target						<< " kg"		<< G4endl
			<< "|   LXe_Gate                         "	<<	Mass_Gate						<< " kg"		<< G4endl
			//<< "|   LXe_Cathode                     "	<<	Mass_Cathode					<< " kg"		<< G4endl
			<< "|   TargetTotal                      "	<<	Mass_TargetTotal		  		<< " kg"		<< G4endl
			<< "| "																							<< G4endl
			<< 	" ============================================================"								<< G4endl;

	// G4cout << G4endl << "The materials defined are : "	<< G4endl << G4endl;	
	//G4cout << *(G4Material::GetMaterialTable())			<< G4endl;


	/////////////////////////
	// Sensitive Detectors //
	///////////////////////////////////////////////////////////////////////////////////////////////////
	G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

	// target volume, xenon sesitivity
	XuerichLXeSensitiveDetector *pLXeSD = new XuerichLXeSensitiveDetector("Xuerich/LXeSD");
	pSDManager->AddNewDetector(pLXeSD);
	LXeTarget_log	->SetSensitiveDetector(pLXeSD);
	LXeGate_log		->SetSensitiveDetector(pLXeSD);
	//LXeCathode_log	->SetSensitiveDetector(pLXeSD);

	// dead xenon
	XuerichDeadXeSensitiveDetector *pDeadXeSD = new XuerichDeadXeSensitiveDetector("Xuerich/DeadXeSD");
	pSDManager->AddNewDetector(pDeadXeSD);
	LXeVol_log->SetSensitiveDetector(pDeadXeSD);

	// gas xenon
	XuerichGXeSensitiveDetector *pGXeSD = new XuerichGXeSensitiveDetector("Xuerich/GXeSD");
	pSDManager->AddNewDetector(pGXeSD);
	GXeVol_log->SetSensitiveDetector(pGXeSD);

	// pmt sensitivity
	XuerichPmtSensitiveDetector *pPmtSD = new XuerichPmtSensitiveDetector("Xuerich/PmtSD");
	pSDManager->AddNewDetector(pPmtSD);
	topPMTcathode_log		->SetSensitiveDetector(pPmtSD);
	bottomPMTcathode_log	->SetSensitiveDetector(pPmtSD);

	// NaI
	XuerichNaISensitiveDetector *pNaISD = new XuerichNaISensitiveDetector("Xuerich/NaISD");
	pSDManager->AddNewDetector(pNaISD);
	NaI_Crystal_log	->SetSensitiveDetector(pNaISD);

	// teflon sensitivity
	XuerichTeflonSensitiveDetector *pTeflonSD = new XuerichTeflonSensitiveDetector("Xuerich/TeflonSD");
	pSDManager->AddNewDetector(pTeflonSD);
	Washer1_log						->SetSensitiveDetector(pTeflonSD);
	Washer2_log						->SetSensitiveDetector(pTeflonSD);
	Washer3_log						->SetSensitiveDetector(pTeflonSD);
	Washer4_log						->SetSensitiveDetector(pTeflonSD);
	Support1_log					->SetSensitiveDetector(pTeflonSD);
	Support2_log					->SetSensitiveDetector(pTeflonSD);
	Support3_log					->SetSensitiveDetector(pTeflonSD);
	Support4_log					->SetSensitiveDetector(pTeflonSD);
	TopClamp_log					->SetSensitiveDetector(pTeflonSD);
	TopHolderUpperTube_log			->SetSensitiveDetector(pTeflonSD);
	TopHolderLowerTube_log			->SetSensitiveDetector(pTeflonSD);
	TopHolderSquare_log				->SetSensitiveDetector(pTeflonSD);
	ExtractionSpacerUpperHalf_log	->SetSensitiveDetector(pTeflonSD);
	ExtractionSpacerLowerHalf_log	->SetSensitiveDetector(pTeflonSD);
	DriftSpacerUpperSquare_log		->SetSensitiveDetector(pTeflonSD);
	DriftSpacerUpperTube_log		->SetSensitiveDetector(pTeflonSD);
	DriftSpacerMiddleTube_log		->SetSensitiveDetector(pTeflonSD);
	DriftSpacerLowerTube_log		->SetSensitiveDetector(pTeflonSD);
	DriftSpacerLowerSquare_log		->SetSensitiveDetector(pTeflonSD);
	BottomHolderUpperTube_log		->SetSensitiveDetector(pTeflonSD);
	BottomHolderMiddleTube_log		->SetSensitiveDetector(pTeflonSD);
	BottomHolderLowerTube_log		->SetSensitiveDetector(pTeflonSD);
	Filler_log						->SetSensitiveDetector(pTeflonSD);
	BottomClampUpperPart_log		->SetSensitiveDetector(pTeflonSD);
	BottomClampMiddlePart_log		->SetSensitiveDetector(pTeflonSD);
	BottomClampLowerPart_log		->SetSensitiveDetector(pTeflonSD);
	topPMTbase_log					->SetSensitiveDetector(pTeflonSD);
	bottomPMTbase_log				->SetSensitiveDetector(pTeflonSD);

	// sensitivity in steel
	XuerichSteelSensitiveDetector *pSteelSD = new XuerichSteelSensitiveDetector("Xuerich/SteelSD");
	pSDManager->AddNewDetector(pSteelSD);
	
	OuterCanTube_log				->SetSensitiveDetector(pSteelSD);
	OuterCanBottom_log				->SetSensitiveDetector(pSteelSD);
	OuterCanLowerFlange_log			->SetSensitiveDetector(pSteelSD);
	OuterCanUpperFlange_log			->SetSensitiveDetector(pSteelSD);
	InnerCanTube_log				->SetSensitiveDetector(pSteelSD);
	InnerCanBottom_log				->SetSensitiveDetector(pSteelSD);
	InnerCanLowerFlange_log			->SetSensitiveDetector(pSteelSD);
	InnerCanUpperFlange_log			->SetSensitiveDetector(pSteelSD);
	SteelHolderUpperSquare_log		->SetSensitiveDetector(pSteelSD);
	SteelHolderLowerSquare_log		->SetSensitiveDetector(pSteelSD);
	SteelHolderTube_log				->SetSensitiveDetector(pSteelSD);
	AnodeGrid_log					->SetSensitiveDetector(pSteelSD);
	CathodeGrid_log					->SetSensitiveDetector(pSteelSD);
	GateGrid_log					->SetSensitiveDetector(pSteelSD);
	topPMTcasing_log				->SetSensitiveDetector(pSteelSD);
	bottomPMTcasing_log				->SetSensitiveDetector(pSteelSD);
	Pipe1in_log						->SetSensitiveDetector(pSteelSD);
	Pipe2in_log						->SetSensitiveDetector(pSteelSD);
	Pipe3in_log						->SetSensitiveDetector(pSteelSD);
	Pipe4in_log						->SetSensitiveDetector(pSteelSD);
	Pipe5in_log						->SetSensitiveDetector(pSteelSD);
	Pipe6in_log						->SetSensitiveDetector(pSteelSD);
	Pipe7in_log						->SetSensitiveDetector(pSteelSD);
	Pipe8in_log						->SetSensitiveDetector(pSteelSD);
	Pipe9in_log						->SetSensitiveDetector(pSteelSD);
	Pipe10in_log					->SetSensitiveDetector(pSteelSD);
	Pipe11in_log					->SetSensitiveDetector(pSteelSD);
	Pipe1out_log					->SetSensitiveDetector(pSteelSD);
	Pipe2out_log					->SetSensitiveDetector(pSteelSD);
	Pipe3out_log					->SetSensitiveDetector(pSteelSD);
	Pipe4out_log					->SetSensitiveDetector(pSteelSD);
	Pipe5out_log					->SetSensitiveDetector(pSteelSD);
	Pipe6out_log					->SetSensitiveDetector(pSteelSD);
	Pipe7out_log					->SetSensitiveDetector(pSteelSD);
	Pipe8out_log					->SetSensitiveDetector(pSteelSD);
	Pipe9out_log					->SetSensitiveDetector(pSteelSD);
	Pipe10out_log					->SetSensitiveDetector(pSteelSD);
	Pipe11out_log					->SetSensitiveDetector(pSteelSD);
	Pipe1flange_log					->SetSensitiveDetector(pSteelSD);
	Pipe2flange_log					->SetSensitiveDetector(pSteelSD);
	Pipe3flange_log					->SetSensitiveDetector(pSteelSD);
	Pipe4flange_log					->SetSensitiveDetector(pSteelSD);
	Pipe5flange1_log				->SetSensitiveDetector(pSteelSD);
	Pipe5flange2_log				->SetSensitiveDetector(pSteelSD);
	Pipe6flange_log					->SetSensitiveDetector(pSteelSD);
	Pipe8flange_log					->SetSensitiveDetector(pSteelSD);
	Pipe9flange_log					->SetSensitiveDetector(pSteelSD);
	Pipe10flange_log				->SetSensitiveDetector(pSteelSD);
	Pipe11flange_log				->SetSensitiveDetector(pSteelSD);

	// sensitivity in aluminum
	XuerichAlSensitiveDetector *pAlSD = new XuerichAlSensitiveDetector("Xuerich/AlSD");
	pSDManager->AddNewDetector(pAlSD);

	XuStand_VertBar_LeftFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_VertBar_LeftBack_log	->SetSensitiveDetector(pAlSD);
	XuStand_VertBar_RightFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_VertBar_RightBack_log	->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_BotLeft_log		->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_BotRight_log	->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_BotFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_BotBack_log		->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_MidLeft_log		->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_MidRight_log	->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_MidFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_MidBack_log		->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_TopLeft_log		->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_TopRight_log	->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_TopFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_HorzBar_TopBack_log		->SetSensitiveDetector(pAlSD);
	XuStand_WingFront_log			->SetSensitiveDetector(pAlSD);
	XuStand_WingBack_log			->SetSensitiveDetector(pAlSD);
	XuStand_BottomPlate_log			->SetSensitiveDetector(pAlSD);
	XuStand_CylinderLeftFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_CylinderLeftBack_log	->SetSensitiveDetector(pAlSD);
	XuStand_CylinderRightFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_CylinderRightBack_log	->SetSensitiveDetector(pAlSD);
	XuStand_CradlePanelLong_log		->SetSensitiveDetector(pAlSD);
	XuStand_CradlePanelShort_log	->SetSensitiveDetector(pAlSD);
	XuStand_CradleBarLeftFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_CradleBarLeftBack_log	->SetSensitiveDetector(pAlSD);
	XuStand_CradleBarRightFront_log	->SetSensitiveDetector(pAlSD);
	XuStand_CradleBarRightBack_log	->SetSensitiveDetector(pAlSD);
	XuStand_Goniometer_log			->SetSensitiveDetector(pAlSD);

	NaIstand_VertBarOut_LeftFront_log	->SetSensitiveDetector(pAlSD);
	NaIstand_VertBarOut_LeftBack_log	->SetSensitiveDetector(pAlSD);
	NaIstand_VertBarOut_RightFront_log	->SetSensitiveDetector(pAlSD);
	NaIstand_VertBarOut_RightBack_log	->SetSensitiveDetector(pAlSD);
	NaIstand_VertBarIn_LeftFront_log	->SetSensitiveDetector(pAlSD);
	NaIstand_VertBarIn_LeftBack_log		->SetSensitiveDetector(pAlSD);
	NaIstand_VertBarIn_RightFront_log	->SetSensitiveDetector(pAlSD);
	NaIstand_VertBarIn_RightBack_log	->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarOut_Left_log		->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarOut_Right_log		->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarOut_Front_log		->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarOut_Back_log		->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarIn_Left_log			->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarIn_Right_log		->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarIn_Front_log		->SetSensitiveDetector(pAlSD);
	NaIstand_HorzBarIn_Back_log			->SetSensitiveDetector(pAlSD);
	NaIstand_PlateTop_log				->SetSensitiveDetector(pAlSD);
	NaIstand_PlateBot_log				->SetSensitiveDetector(pAlSD);
	NaIstand_CradleBarLeft_log			->SetSensitiveDetector(pAlSD);
	NaIstand_CradleBarRight_log			->SetSensitiveDetector(pAlSD);
	NaIstand_CradlePanel_log			->SetSensitiveDetector(pAlSD);
	NaI_CrystalHousingTop_log			->SetSensitiveDetector(pAlSD);
	NaI_CrystalHousingBot_log			->SetSensitiveDetector(pAlSD);

	LeadChannelStand_VertBarOut_LeftFront_log	->SetSensitiveDetector(pAlSD);
	LeadChannelStand_VertBarOut_LeftBack_log	->SetSensitiveDetector(pAlSD);
	LeadChannelStand_VertBarOut_RightFront_log	->SetSensitiveDetector(pAlSD);
	LeadChannelStand_VertBarOut_RightBack_log	->SetSensitiveDetector(pAlSD);
	LeadChannelStand_VertBarIn_LeftFront_log	->SetSensitiveDetector(pAlSD);
	LeadChannelStand_VertBarIn_LeftBack_log		->SetSensitiveDetector(pAlSD);
	LeadChannelStand_VertBarIn_RightFront_log	->SetSensitiveDetector(pAlSD);
	LeadChannelStand_VertBarIn_RightBack_log	->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarOut_Left_log		->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarOut_Right_log		->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarOut_Front_log		->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarOut_Back_log		->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarIn_Left_log			->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarIn_Right_log		->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarIn_Front_log		->SetSensitiveDetector(pAlSD);
	LeadChannelStand_HorzBarIn_Back_log			->SetSensitiveDetector(pAlSD);
	LeadChannelStand_PlateTop_log				->SetSensitiveDetector(pAlSD);
	LeadChannelStand_PlateBot_log				->SetSensitiveDetector(pAlSD);

	ApertureFrame_log					->SetSensitiveDetector(pAlSD);
	
	AlCanTube_log						->SetSensitiveDetector(pAlSD);
	AlCanBottom_log						->SetSensitiveDetector(pAlSD);
	
	// sensitivity in LEAD
	XuerichLeadSensitiveDetector *pLeadSD = new XuerichLeadSensitiveDetector("Xuerich/LeadSD");
	pSDManager->AddNewDetector(pLeadSD);
	NaIcastle_log			->SetSensitiveDetector(pLeadSD);
	Collimator_log			->SetSensitiveDetector(pLeadSD);
	LeadChannel_log			->SetSensitiveDetector(pLeadSD);
	Aperture_log			->SetSensitiveDetector(pLeadSD);
	
	// sensitivity IN THE REST OF MATERIALS
	XuerichRestSensitiveDetector *pRestSD = new XuerichRestSensitiveDetector("Xuerich/RestSD");
	pSDManager->AddNewDetector(pRestSD);
	topPMTwindow_log		->SetSensitiveDetector(pRestSD);
	bottomPMTwindow_log		->SetSensitiveDetector(pRestSD);
	topPMTcathode_log		->SetSensitiveDetector(pRestSD);
	bottomPMTcathode_log	->SetSensitiveDetector(pRestSD);
	CopperFinger_log		->SetSensitiveDetector(pRestSD);
	NaI_LightShieldTop_log	->SetSensitiveDetector(pRestSD);
	NaI_LightShieldBot_log	->SetSensitiveDetector(pRestSD);
	NaI_LightShieldMid_log	->SetSensitiveDetector(pRestSD);
	DewarOutTube_log		->SetSensitiveDetector(pRestSD);
	DewarOutBot_log			->SetSensitiveDetector(pRestSD);
	DewarInTube_log			->SetSensitiveDetector(pRestSD);
	DewarInBot_log			->SetSensitiveDetector(pRestSD);
	DewarTop_log			->SetSensitiveDetector(pRestSD);

			


//--------------	END PROGRAM	------------------------
return Laboratory_phys;
//------------------------------------------------------

}



  ////////////////////////
 // Optical parameters //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void XuerichDetectorConstruction::SetTeflonReflectivity(G4double dReflectivity)
{
	//G4Material *Teflon = G4Material::GetMaterial(G4String("Teflon"));

	if(Teflon)
	{
		G4cout << "\n----> Setting Teflon reflectivity to " << dReflectivity << G4endl;

		G4MaterialPropertiesTable *pTeflonPropertiesTable = Teflon->GetMaterialPropertiesTable();
		
		G4double teflon_PP[] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
		G4double teflon_REFL[] = {dReflectivity, dReflectivity, dReflectivity};
		pTeflonPropertiesTable->RemoveProperty("REFLECTIVITY");
		pTeflonPropertiesTable->AddProperty("REFLECTIVITY", teflon_PP, teflon_REFL, 3);
	}
	else
	{
		G4cout << "!!!!> Teflon material not found!" << G4endl;
		exit(-1);
	}
}

void XuerichDetectorConstruction::SetSS304LSteelReflectivity(G4double dSteelReflectivity)
{
	//G4Material *SSteel = G4Material::GetMaterial(G4String("SS304LSteel"));

	if(SSteel)
	{
		G4cout << "\n----> Setting SS304LSteel reflectivity to " << dSteelReflectivity << G4endl;

		G4MaterialPropertiesTable *pSteelPropertiesTable = SSteel->GetMaterialPropertiesTable();
		
		G4double Steel_PP[] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
		G4double Steel_REFL[] = {dSteelReflectivity, dSteelReflectivity, dSteelReflectivity};
		pSteelPropertiesTable->RemoveProperty("REFLECTIVITY");
		pSteelPropertiesTable->AddProperty("REFLECTIVITY", Steel_PP, Steel_REFL, 3);
	}
	else
	{
		G4cout << "!!!!> SS304LSteel material not found!" << G4endl;
		exit(-1);
	}
}


void XuerichDetectorConstruction::SetLXeScintillation(G4bool bScintillation)
{
	G4cout << "----> Setting LXe(GXe) scintillation to " << bScintillation << G4endl;
			
	G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));
	if(pLXeMaterial)
	{	
		G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
		if(bScintillation)
			pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1000./(1.0*keV));
	}
	else
	{
		G4cout << "ls!> LXe materials not found!" << G4endl;
		exit(-1);
	}
	
	G4Material *pGXeMaterial = G4Material::GetMaterial(G4String("GXe"));
	if(pGXeMaterial)
	{	
		G4MaterialPropertiesTable *pGXePropertiesTable = pGXeMaterial->GetMaterialPropertiesTable();
		if(bScintillation)
			pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1000/(1.0*keV));
	}
	else
	{
		G4cout << "ls!> GXe materials not found!" << G4endl;
		exit(-1);
	}


}



void XuerichDetectorConstruction::SetLXeAbsorbtionLength(G4double dAbsorbtionLength)
{
	G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));

	if(pLXeMaterial)
	{
		G4cout << "----> Setting LXe absorbtion length to " << dAbsorbtionLength/cm << " cm" << G4endl;

		G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
			
			G4double LXe_PP[] = {6.91*eV, 6.98*eV, 7.05*eV};
			G4double LXe_ABSL[] = {dAbsorbtionLength, dAbsorbtionLength, dAbsorbtionLength};
			pLXePropertiesTable->RemoveProperty("ABSLENGTH");
			pLXePropertiesTable->AddProperty("ABSLENGTH", LXe_PP, LXe_ABSL, 3);
	}
	else
	{
		G4cout << "ls!> LXe materials not found!" << G4endl;
		exit(-1);
	}
}

void XuerichDetectorConstruction::SetLXeRayScatterLength(G4double dRayScatterLength)
{
  G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));
  
  if(pLXeMaterial)
    {
      
      G4cout << "----> Setting LXe scattering length to " << dRayScatterLength/cm << " cm" << G4endl;
      
      G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
              
               G4double LXe_PP[] = {6.91*eV, 6.98*eV, 7.05*eV};
               G4double LXe_SCAT[] = {dRayScatterLength, dRayScatterLength, dRayScatterLength};
               pLXePropertiesTable->RemoveProperty("RAYLEIGH");
               pLXePropertiesTable->AddProperty("RAYLEIGH", LXe_PP, LXe_SCAT, 3); 
    }
  else
    {
      G4cout << "ls!> LXe materials not found!" << G4endl;
      exit(-1);
    }

}




