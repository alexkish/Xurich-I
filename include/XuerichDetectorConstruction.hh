/////////////////////////////////////////////////////////
//													   //
//   ++++++++++++++++++++++++++++++++++++++++++++++++  //
//   + Alexander Kish for XENON-XUERICH experiment	+  //
//   + UZH, 2008									+  //
// 	 ++++++++++++++++++++++++++++++++++++++++++++++++  //
//													   //
/////////////////////////////////////////////////////////
#ifndef XuerichDetectorConstruction_h
#define XuerichDetectorConstruction_h 1

#include <globals.hh>

#include <vector>
#include <map>


#include "G4VUserDetectorConstruction.hh"
#include "XuerichDetectorConstruction.hh"


class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VisAttributes;
class G4Material;
class G4MaterialPropertiesTable;
class G4SubtractionSolid;
class G4UnionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

class XuerichDetectorMessenger;


class XuerichDetectorConstruction: public G4VUserDetectorConstruction
{

  public:
    XuerichDetectorConstruction();
   ~XuerichDetectorConstruction();

    G4VPhysicalVolume *Construct();

	void SetTeflonReflectivity(G4double dReflectivity);
	void SetSS304LSteelReflectivity(G4double dReflectivity);
	void SetLXeScintillation(G4bool dScintillation);
	void SetLXeAbsorbtionLength(G4double dAbsorbtionLength);
    void SetLXeRayScatterLength(G4double dRayScatterLength);

	   
  private:

	G4Material *	Air;
	G4Material *	Vacuum;
	G4Material *	SSteel;
	G4Material *	MuMetal;
	G4Material *	Teflon;
	G4Material *	LXe;
	G4Material *	GXe;
	//G4Material *	LiqXe;
	//G4Material *	GasXe;
	G4Material *	LAr;
	G4Material *	LiqAr;
	G4Material *	GasAr;
	G4Material *	Quartz;
	G4Material *	metalAl;
	G4Material *	Alwire;
	G4Material *	PhotoCathodeAluminium;
	G4Material *	Lead;
	G4Material *	Poly;
	G4Material *	Copper;
	G4Material *	Concrete;
	G4Material *	BPE;
	G4Material *	Paraffin;
	G4Material *	Lead_m;
	G4Material *	NaI;

	//G4MaterialPropertiesTable*      teflon_prop;




	////////////////////////////////////////////////////
	// Logical Volumes
	////////////////////////////////////////////////////
	G4LogicalVolume *	Laboratory_log;	
	// xuerich stand
	G4LogicalVolume	*	XuStand_VertBar_LeftFront_log;
	G4LogicalVolume	*	XuStand_VertBar_LeftBack_log;
	G4LogicalVolume	*	XuStand_VertBar_RightFront_log;
	G4LogicalVolume	*	XuStand_VertBar_RightBack_log;

	G4LogicalVolume	*	XuStand_HorzBar_BotLeft_log;
	G4LogicalVolume	*	XuStand_HorzBar_BotRight_log;
	G4LogicalVolume	*	XuStand_HorzBar_BotFront_log;
	G4LogicalVolume	*	XuStand_HorzBar_BotBack_log;

	G4LogicalVolume	*	XuStand_HorzBar_MidLeft_log;
	G4LogicalVolume	*	XuStand_HorzBar_MidRight_log;
	G4LogicalVolume	*	XuStand_HorzBar_MidFront_log;
	G4LogicalVolume	*	XuStand_HorzBar_MidBack_log;

	G4LogicalVolume	*	XuStand_HorzBar_TopLeft_log;
	G4LogicalVolume	*	XuStand_HorzBar_TopRight_log;
	G4LogicalVolume	*	XuStand_HorzBar_TopFront_log;
	G4LogicalVolume	*	XuStand_HorzBar_TopBack_log;

	G4LogicalVolume	*	XuStand_WingFront_log;
	G4LogicalVolume	*	XuStand_WingBack_log;

	G4LogicalVolume	*	XuStand_BottomPlate_log;

	G4LogicalVolume	* 	XuStand_CylinderLeftFront_log;	
	G4LogicalVolume	* 	XuStand_CylinderLeftBack_log;	
	G4LogicalVolume	* 	XuStand_CylinderRightFront_log;	
	G4LogicalVolume	* 	XuStand_CylinderRightBack_log;	

	G4LogicalVolume	*  	XuStand_CradlePanelLong_log;
	G4LogicalVolume	*	XuStand_CradlePanelShort_log;

	G4LogicalVolume	*	XuStand_CradleBarLeftFront_log;
	G4LogicalVolume	*	XuStand_CradleBarLeftBack_log;
	G4LogicalVolume	*	XuStand_CradleBarRightFront_log;
	G4LogicalVolume	*	XuStand_CradleBarRightBack_log;

	G4LogicalVolume	*	XuStand_Goniometer_log;
	// Aluminum can
	G4LogicalVolume	*	AlCanTube_log;
	G4LogicalVolume	*	AlCanBottom_log;
	// copper finger
	G4LogicalVolume	*	CopperFinger_log;
	// dewar
	G4LogicalVolume	*	DewarOutTube_log;
	G4LogicalVolume	*	DewarOutBot_log;
	G4LogicalVolume	*	DewarInTube_log;
	G4LogicalVolume	*	DewarInBot_log;
	G4LogicalVolume	*	DewarTop_log;
	// cryostat
	G4LogicalVolume *	OuterCanTube_log;
	G4LogicalVolume *	OuterCanBottom_log;
	G4LogicalVolume *	OuterCanLowerFlange_log;
	G4LogicalVolume *	OuterCanUpperFlange_log;

	G4LogicalVolume *	Pipe1in_log;
	G4LogicalVolume *	Pipe2in_log;
	G4LogicalVolume *	Pipe3in_log;
	G4LogicalVolume *	Pipe4in_log;
	G4LogicalVolume *	Pipe5in_log;
	G4LogicalVolume *	Pipe6in_log;
	G4LogicalVolume *	Pipe7in_log;
	G4LogicalVolume *	Pipe8in_log;
	G4LogicalVolume *	Pipe9in_log;
	G4LogicalVolume *	Pipe10in_log;
	G4LogicalVolume *	Pipe11in_log;

	G4LogicalVolume *	Pipe1out_log;
	G4LogicalVolume *	Pipe2out_log;
	G4LogicalVolume *	Pipe3out_log;
	G4LogicalVolume *	Pipe4out_log;
	G4LogicalVolume *	Pipe5out_log;
	G4LogicalVolume *	Pipe6out_log;
	G4LogicalVolume *	Pipe7out_log;
	G4LogicalVolume *	Pipe8out_log;
	G4LogicalVolume *	Pipe9out_log;
	G4LogicalVolume *	Pipe10out_log;
	G4LogicalVolume *	Pipe11out_log;

	G4LogicalVolume *	Pipe1flange_log;
	G4LogicalVolume *	Pipe2flange_log;
	G4LogicalVolume *	Pipe3flange_log;
	G4LogicalVolume *	Pipe4flange_log;
	G4LogicalVolume *	Pipe5flange1_log;
	G4LogicalVolume *	Pipe5flange2_log;
	G4LogicalVolume *	Pipe6flange_log;
	G4LogicalVolume *	Pipe8flange_log;
	G4LogicalVolume *	Pipe9flange_log;
	G4LogicalVolume *	Pipe10flange_log;
	G4LogicalVolume *	Pipe11flange_log;

	G4LogicalVolume *	InnerCanTube_log;
	G4LogicalVolume *	InnerCanBottom_log;
	G4LogicalVolume *	InnerCanLowerFlange_log;
	G4LogicalVolume *	InnerCanUpperFlange_log;

	G4LogicalVolume *	SteelHolderTube_log;	
	G4LogicalVolume *	SteelHolderUpperSquare_log;
	G4LogicalVolume *	SteelHolderLowerSquare_log;

	G4LogicalVolume *	Support1_log;	////////////											
	G4LogicalVolume *	Support2_log;	//	1 - 2 //
	G4LogicalVolume *	Support3_log;	//	|	| //
	G4LogicalVolume *	Support4_log;	//	3 - 4 //
										////////////
	G4LogicalVolume *	Washer1_log;	////////////
	G4LogicalVolume *	Washer2_log;	//	1 - 2 //
	G4LogicalVolume *	Washer3_log;	//	|	| //
	G4LogicalVolume *	Washer4_log;	//	3 - 4 //
										////////////
	G4LogicalVolume *	TopClamp_log;
	
	G4LogicalVolume	*	TopHolderUpperTube_log;	
	G4LogicalVolume	*	TopHolderLowerTube_log;	
	G4LogicalVolume	*	TopHolderSquare_log;	

	G4LogicalVolume *	ExtractionSpacerUpperHalf_log;
	G4LogicalVolume *	ExtractionSpacerLowerHalf_log;
	
	G4LogicalVolume *	DriftSpacerUpperSquare_log;
	G4LogicalVolume *	DriftSpacerUpperTube_log;
	G4LogicalVolume *	DriftSpacerMiddleTube_log;
	G4LogicalVolume *	DriftSpacerLowerTube_log;
	G4LogicalVolume *	DriftSpacerLowerSquare_log;
	
	G4LogicalVolume	*	BottomHolderUpperTube_log;
	G4LogicalVolume	*	BottomHolderMiddleTube_log;
	G4LogicalVolume	*	BottomHolderLowerTube_log;
	
	G4LogicalVolume *	BottomClampUpperPart_log;
	G4LogicalVolume *	BottomClampMiddlePart_log;
	G4LogicalVolume *	BottomClampLowerPart_log;

	G4LogicalVolume	*	Filler_log;
	
	G4LogicalVolume *	topPMTbase_log;
	G4LogicalVolume *	topPMTcasing_log;
	G4LogicalVolume *	topPMTinterior_log;
	G4LogicalVolume *	topPMTwindow_log;
	G4LogicalVolume *	topPMTcathode_log;
	
	G4LogicalVolume *	bottomPMTbase_log;
	G4LogicalVolume *	bottomPMTcasing_log;
	G4LogicalVolume *	bottomPMTinterior_log;
	G4LogicalVolume *	bottomPMTwindow_log;
	G4LogicalVolume *	bottomPMTcathode_log;
	
	G4LogicalVolume *	AnodeGrid_log;
	G4LogicalVolume *	CathodeGrid_log;
	G4LogicalVolume *	GateGrid_log;

	G4LogicalVolume	*	LXeVol_log;
	G4LogicalVolume	*	GXeVol_log;
	G4LogicalVolume	*	LXeTarget_log;
	G4LogicalVolume	*	LXeGate_log;
	G4LogicalVolume	*	LXeCathode_log;
	// NaI
	G4LogicalVolume	*NaI_Crystal_log;
	G4LogicalVolume	*NaI_CrystalHousingTop_log;
	G4LogicalVolume	*NaI_CrystalHousingBot_log;
	G4LogicalVolume	*NaI_LightShieldTop_log;
	G4LogicalVolume	*NaI_LightShieldBot_log;
	G4LogicalVolume	*NaI_LightShieldMid_log;
	//G4LogicalVolume	*NaI_LightShieldVac_log;
	//G4LogicalVolume	*NaI_PMTcasing_log;

	G4LogicalVolume	*NaIstand_VertBarOut_LeftFront_log;
	G4LogicalVolume	*NaIstand_VertBarOut_LeftBack_log;	
	G4LogicalVolume	*NaIstand_VertBarOut_RightFront_log;	
	G4LogicalVolume	*NaIstand_VertBarOut_RightBack_log;	
	
	G4LogicalVolume	*NaIstand_VertBarIn_LeftFront_log;	
	G4LogicalVolume	*NaIstand_VertBarIn_LeftBack_log;	
	G4LogicalVolume	*NaIstand_VertBarIn_RightFront_log;	
	G4LogicalVolume	*NaIstand_VertBarIn_RightBack_log;	

	G4LogicalVolume	*NaIstand_HorzBarOut_Left_log;		
	G4LogicalVolume	*NaIstand_HorzBarOut_Right_log;		
	G4LogicalVolume	*NaIstand_HorzBarOut_Front_log;		
	G4LogicalVolume	*NaIstand_HorzBarOut_Back_log;		

	G4LogicalVolume	*NaIstand_HorzBarIn_Left_log;		
	G4LogicalVolume	*NaIstand_HorzBarIn_Right_log;		
	G4LogicalVolume	*NaIstand_HorzBarIn_Front_log;		
	G4LogicalVolume	*NaIstand_HorzBarIn_Back_log;		

	G4LogicalVolume	*NaIstand_PlateTop_log;				
	G4LogicalVolume	*NaIstand_PlateBot_log;				

	G4LogicalVolume	*NaIstand_CradleBarLeft_log;				
	G4LogicalVolume	*NaIstand_CradleBarRight_log;				
	G4LogicalVolume	*NaIstand_CradlePanel_log;				

	G4LogicalVolume	*LeadChannelStand_VertBarOut_LeftFront_log;
	G4LogicalVolume	*LeadChannelStand_VertBarOut_LeftBack_log;	
	G4LogicalVolume	*LeadChannelStand_VertBarOut_RightFront_log;	
	G4LogicalVolume	*LeadChannelStand_VertBarOut_RightBack_log;	
	
	G4LogicalVolume	*LeadChannelStand_VertBarIn_LeftFront_log;	
	G4LogicalVolume	*LeadChannelStand_VertBarIn_LeftBack_log;	
	G4LogicalVolume	*LeadChannelStand_VertBarIn_RightFront_log;	
	G4LogicalVolume	*LeadChannelStand_VertBarIn_RightBack_log;	

	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Left_log;		
	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Right_log;		
	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Front_log;		
	G4LogicalVolume	*LeadChannelStand_HorzBarOut_Back_log;		

	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Left_log;		
	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Right_log;		
	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Front_log;		
	G4LogicalVolume	*LeadChannelStand_HorzBarIn_Back_log;		

	G4LogicalVolume	*LeadChannelStand_PlateTop_log;				
	G4LogicalVolume	*LeadChannelStand_PlateBot_log;				

	G4LogicalVolume	*NaIcastle_log;
	G4LogicalVolume	*Aperture_log;
	G4LogicalVolume	*ApertureFrame_log;

	G4LogicalVolume	*Collimator_log;

	G4LogicalVolume	*LeadChannel_log;
	
	////////////////////////////////////////////////////////////////////
	// PHYSICAL VOLUMES
	////////////////////////////////////////////////////////////////////
	G4VPhysicalVolume *	Laboratory_phys;

	G4VPhysicalVolume *	XuStand_VertBar_LeftFront_phys;
	G4VPhysicalVolume *	XuStand_VertBar_LeftBack_phys;
	G4VPhysicalVolume *	XuStand_VertBar_RightFront_phys;
	G4VPhysicalVolume *	XuStand_VertBar_RightBack_phys;

	G4VPhysicalVolume *	XuStand_HorzBar_BotLeft_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_BotRight_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_BotFront_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_BotBack_phys;

	G4VPhysicalVolume *	XuStand_HorzBar_MidLeft_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_MidRight_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_MidFront_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_MidBack_phys;

	G4VPhysicalVolume *	XuStand_HorzBar_TopLeft_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_TopRight_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_TopFront_phys;
	G4VPhysicalVolume *	XuStand_HorzBar_TopBack_phys;

	G4VPhysicalVolume *	XuStand_BottomPlate_phys;

	G4VPhysicalVolume *	XuStand_WingFront_phys;
	G4VPhysicalVolume *	XuStand_WingBack_phys;

	G4VPhysicalVolume * XuStand_CylinderLeftFront_phys;	
	G4VPhysicalVolume * XuStand_CylinderLeftBack_phys;	
	G4VPhysicalVolume * XuStand_CylinderRightFront_phys;	
	G4VPhysicalVolume * XuStand_CylinderRightBack_phys;

	G4VPhysicalVolume * XuStand_CradlePanelLong_phys;
	G4VPhysicalVolume *	XuStand_CradlePanelShort_phys;

	G4VPhysicalVolume *	XuStand_CradleBarLeftFront_phys;
	G4VPhysicalVolume *	XuStand_CradleBarLeftBack_phys;
	G4VPhysicalVolume *	XuStand_CradleBarRightFront_phys;
	G4VPhysicalVolume *	XuStand_CradleBarRightBack_phys;

	G4VPhysicalVolume *	XuStand_Goniometer_phys;

	G4VPhysicalVolume *	AlCanTube_phys;
	G4VPhysicalVolume *	AlCanBotom_phys;
	G4VPhysicalVolume *	CopperFinger_phys;

	G4VPhysicalVolume *	DewarOutTube_phys;
	G4VPhysicalVolume *	DewarOutBot_phys;
	G4VPhysicalVolume *	DewarInTube_phys;
	G4VPhysicalVolume *	DewarInBot_phys;
	G4VPhysicalVolume *	DewarTop_phys;
	
	G4VPhysicalVolume *	OuterCanTube_phys;
	G4VPhysicalVolume *	OuterCanBottom_phys;
	G4VPhysicalVolume *	OuterCanLowerFlange_phys;
	G4VPhysicalVolume *	OuterCanUpperFlange_phys;

	G4VPhysicalVolume *	Pipe1in_phys;
	G4VPhysicalVolume *	Pipe2in_phys;
	G4VPhysicalVolume *	Pipe3in_phys;
	G4VPhysicalVolume *	Pipe4in_phys;
	G4VPhysicalVolume *	Pipe5in_phys;
	G4VPhysicalVolume *	Pipe6in_phys;
	G4VPhysicalVolume *	Pipe7in_phys;
	G4VPhysicalVolume *	Pipe8in_phys;
	G4VPhysicalVolume *	Pipe9in_phys;
	G4VPhysicalVolume *	Pipe10in_phys;
	G4VPhysicalVolume *	Pipe11in_phys;

	G4VPhysicalVolume *	Pipe1out_phys;
	G4VPhysicalVolume *	Pipe2out_phys;
	G4VPhysicalVolume *	Pipe3out_phys;
	G4VPhysicalVolume *	Pipe4out_phys;
	G4VPhysicalVolume *	Pipe5out_phys;
	G4VPhysicalVolume *	Pipe6out_phys;
	G4VPhysicalVolume *	Pipe7out_phys;
	G4VPhysicalVolume *	Pipe8out_phys;
	G4VPhysicalVolume *	Pipe9out_phys;
	G4VPhysicalVolume *	Pipe10out_phys;
	G4VPhysicalVolume *	Pipe11out_phys;

	G4VPhysicalVolume *	Pipe1flange_phys;
	G4VPhysicalVolume *	Pipe2flange_phys;
	G4VPhysicalVolume *	Pipe3flange_phys;
	G4VPhysicalVolume *	Pipe4flange_phys;
	G4VPhysicalVolume *	Pipe5flange1_phys;
	G4VPhysicalVolume *	Pipe5flange2_phys;
	G4VPhysicalVolume *	Pipe6flange_phys;
	G4VPhysicalVolume *	Pipe8flange_phys;
	G4VPhysicalVolume *	Pipe9flange_phys;
	G4VPhysicalVolume *	Pipe10flange_phys;
	G4VPhysicalVolume *	Pipe11flange_phys;

	G4VPhysicalVolume *	InnerCanTube_phys;
	G4VPhysicalVolume *	InnerCanBottom_phys;
	G4VPhysicalVolume *	InnerCanLowerFlange_phys;
	G4VPhysicalVolume *	InnerCanUpperFlange_phys;
	
	G4VPhysicalVolume *	SteelHolderTube_phys;
	G4VPhysicalVolume *	SteelHolderUpperSquare_phys;
	G4VPhysicalVolume *	SteelHolderLowerSquare_phys;
	
	G4VPhysicalVolume *	Support1_phys;	////////////											
	G4VPhysicalVolume *	Support2_phys;	//	1 - 2 //
	G4VPhysicalVolume *	Support3_phys;	//	|	| //
	G4VPhysicalVolume *	Support4_phys;	//	3 - 4 //
										////////////
	G4VPhysicalVolume *	Washer1_phys;	////////////
	G4VPhysicalVolume *	Washer2_phys;	//	1 - 2 //
	G4VPhysicalVolume *	Washer3_phys;	//	|	| //
	G4VPhysicalVolume *	Washer4_phys;	//	3 - 4 //
										////////////
	G4VPhysicalVolume *	TopClamp_phys;
	
	G4VPhysicalVolume *	TopHolderUpperTube_phys;	
	G4VPhysicalVolume *	TopHolderLowerTube_phys;	
	G4VPhysicalVolume *	TopHolderSquare_phys;	
	
	G4VPhysicalVolume *	ExtractionSpacerUpperHalf_phys;
	G4VPhysicalVolume *	ExtractionSpacerLowerHalf_phys;
	
	G4VPhysicalVolume *	DriftSpacerUpperSquare_phys;
	G4VPhysicalVolume *	DriftSpacerUpperTube_phys;
	G4VPhysicalVolume *	DriftSpacerMiddleTube_phys;
	G4VPhysicalVolume *	DriftSpacerLowerTube_phys;
	G4VPhysicalVolume *	DriftSpacerLowerSquare_phys;
	
	G4VPhysicalVolume *	BottomHolderUpperTube_phys;
	G4VPhysicalVolume *	BottomHolderMiddleTube_phys;
	G4VPhysicalVolume *	BottomHolderLowerTube_phys;

	G4VPhysicalVolume *	BottomClampUpperPart_phys;
	G4VPhysicalVolume *	BottomClampMiddlePart_phys;
	G4VPhysicalVolume *	BottomClampLowerPart_phys;

	G4VPhysicalVolume *	Filler_phys;
	
	G4VPhysicalVolume *	topPMTbase_phys;
	G4VPhysicalVolume *	topPMTcasing_phys;
	G4VPhysicalVolume *	topPMTinterior_phys;
	G4VPhysicalVolume *	topPMTwindow_phys;
	G4VPhysicalVolume *	topPMTcathode_phys;
	
	G4VPhysicalVolume *	bottomPMTbase_phys;
	G4VPhysicalVolume *	bottomPMTcasing_phys;
	G4VPhysicalVolume *	bottomPMTinterior_phys;
	G4VPhysicalVolume *	bottomPMTwindow_phys;
	G4VPhysicalVolume *	bottomPMTcathode_phys;
	
	G4VPhysicalVolume *	AnodeGrid_phys;
	G4VPhysicalVolume *	CathodeGrid_phys;
	G4VPhysicalVolume *	GateGrid_phys;

	G4VPhysicalVolume *	LXeVol_phys;
	G4VPhysicalVolume *	GXeVol_phys;
	G4VPhysicalVolume *	LXeTarget_phys;
	G4VPhysicalVolume *	LXeGate_phys;
	G4VPhysicalVolume *	LXeCathode_phys;

	G4VPhysicalVolume * NaI_Crystal_phys;
	G4VPhysicalVolume * NaI_CrystalHousingTop_phys;
	G4VPhysicalVolume * NaI_CrystalHousingBot_phys;
	G4VPhysicalVolume * NaI_LightShieldTop_phys;
	G4VPhysicalVolume * NaI_LightShieldBot_phys;
	G4VPhysicalVolume * NaI_LightShieldMid_phys;
	//G4VPhysicalVolume * NaI_LightShieldVac_phys;
	//G4VPhysicalVolume * NaI_PMTcasing_phys;

	G4VPhysicalVolume	*NaIstand_VertBarOut_LeftFront_phys;	
	G4VPhysicalVolume	*NaIstand_VertBarOut_LeftBack_phys;	
	G4VPhysicalVolume	*NaIstand_VertBarOut_RightFront_phys;	
	G4VPhysicalVolume	*NaIstand_VertBarOut_RightBack_phys;	
	
	G4VPhysicalVolume	*NaIstand_VertBarIn_LeftFront_phys;	
	G4VPhysicalVolume	*NaIstand_VertBarIn_LeftBack_phys;	
	G4VPhysicalVolume	*NaIstand_VertBarIn_RightFront_phys;	
	G4VPhysicalVolume	*NaIstand_VertBarIn_RightBack_phys;	

	G4VPhysicalVolume	*NaIstand_HorzBarOut_Left_phys;		
	G4VPhysicalVolume	*NaIstand_HorzBarOut_Right_phys;		
	G4VPhysicalVolume	*NaIstand_HorzBarOut_Front_phys;		
	G4VPhysicalVolume	*NaIstand_HorzBarOut_Back_phys;		

	G4VPhysicalVolume	*NaIstand_HorzBarIn_Left_phys;		
	G4VPhysicalVolume	*NaIstand_HorzBarIn_Right_phys;		
	G4VPhysicalVolume	*NaIstand_HorzBarIn_Front_phys;		
	G4VPhysicalVolume	*NaIstand_HorzBarIn_Back_phys;		

	G4VPhysicalVolume	*NaIstand_PlateTop_phys;				
	G4VPhysicalVolume	*NaIstand_PlateBot_phys;				

	G4VPhysicalVolume	*NaIstand_CradleBarLeft_phys;				
	G4VPhysicalVolume	*NaIstand_CradleBarRight_phys;				
	G4VPhysicalVolume	*NaIstand_CradlePanel_phys;				

	G4VPhysicalVolume	*LeadChannelStand_VertBarOut_LeftFront_phys;	
	G4VPhysicalVolume	*LeadChannelStand_VertBarOut_LeftBack_phys;	
	G4VPhysicalVolume	*LeadChannelStand_VertBarOut_RightFront_phys;	
	G4VPhysicalVolume	*LeadChannelStand_VertBarOut_RightBack_phys;	
	
	G4VPhysicalVolume	*LeadChannelStand_VertBarIn_LeftFront_phys;	
	G4VPhysicalVolume	*LeadChannelStand_VertBarIn_LeftBack_phys;	
	G4VPhysicalVolume	*LeadChannelStand_VertBarIn_RightFront_phys;	
	G4VPhysicalVolume	*LeadChannelStand_VertBarIn_RightBack_phys;	

	G4VPhysicalVolume	*LeadChannelStand_HorzBarOut_Left_phys;		
	G4VPhysicalVolume	*LeadChannelStand_HorzBarOut_Right_phys;		
	G4VPhysicalVolume	*LeadChannelStand_HorzBarOut_Front_phys;		
	G4VPhysicalVolume	*LeadChannelStand_HorzBarOut_Back_phys;		

	G4VPhysicalVolume	*LeadChannelStand_HorzBarIn_Left_phys;		
	G4VPhysicalVolume	*LeadChannelStand_HorzBarIn_Right_phys;		
	G4VPhysicalVolume	*LeadChannelStand_HorzBarIn_Front_phys;		
	G4VPhysicalVolume	*LeadChannelStand_HorzBarIn_Back_phys;		

	G4VPhysicalVolume	*LeadChannelStand_PlateTop_phys;				
	G4VPhysicalVolume	*LeadChannelStand_PlateBot_phys;				

	G4VPhysicalVolume	*NaIcastle_phys;
	G4VPhysicalVolume	*Aperture_phys;
	G4VPhysicalVolume	*ApertureFrame_phys;

	G4VPhysicalVolume	*Collimator_phys;				

	G4VPhysicalVolume	*LeadChannel_phys;				



	////////////////////////////////////////////////////////////////
	// Visual Attributes
	////////////////////////////////////////////////////////////////
	G4VisAttributes *	Steel_vis;
	G4VisAttributes *	Al_vis;
	G4VisAttributes *	Cu_vis;
	G4VisAttributes *	Teflon_vis;
	G4VisAttributes *	Cryostat_vis;
	G4VisAttributes *	SteelHolder_vis;
	G4VisAttributes *	Grids_vis;
	G4VisAttributes *	PMTcasing_vis;
	G4VisAttributes *	PMTinterior_vis;
	G4VisAttributes *	PMTwindow_vis;
	G4VisAttributes *	PMTcathode_vis;

	G4VisAttributes *	NaI_Crystal_vis;
	G4VisAttributes *	NaI_CrystalHousing_vis;
	G4VisAttributes *	NaI_LightShield_vis;
	G4VisAttributes *	NaI_PMTcasing_vis;

	G4VisAttributes *	NaIcastle_vis;

	G4VisAttributes *	Collimator_vis;

	G4VisAttributes *	LeadChannel_vis;

	
	XuerichDetectorMessenger *m_pDetectorMessenger;

};

#endif

