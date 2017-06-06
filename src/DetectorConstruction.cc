// $Id: DetectorConstruction.cc 2017-05-29 vega $
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"


//G4ThreadLocal

// The Constructor
DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(),
					       fAbsorberPV(0),
					       fTargetPV(0),
					       fCheckOverlaps(true)
{
  // Define Materials
  DefineMaterials();

  // Define parameters
  ComputeParams();

}

// The Destructor
DetectorConstruction::~DetectorConstruction()
{
}


// Now define Material in a function
void DetectorConstruction::DefineMaterials()
{
  G4NistManager* nistManager = G4NistManager::Instance();

  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 

  nistManager->FindOrBuildMaterial("G4_Tc");
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  nistManager->FindOrBuildMaterial("G4_AIR");
  nistManager->FindOrBuildMaterial("G4_Galactic");

  // Print the material synthesized
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}


// Compute Parameters
void DetectorConstruction::ComputeParams(){

  // Geometry operations  
  RotXP->rotateX( 90*deg);   RotYP->rotateY( 90*deg);  RotZP->rotateZ( 90*deg);
  RotXM->rotateX(-90*deg);   RotYM->rotateY(-90*deg);  RotZM->rotateZ(-90*deg);

  // Geometry Parameters
  // brick full length
  brickX      = 20.0*cm;               brickY        = 10.0*cm;           brickZ       = 5.0*cm;
  HalfBrickX  = 0.5 * brickX;          HalfBrickY    = 0.5 * brickY;      HalfBrickZ   = 0.5 * brickZ;
  b2bGap      = 0.1*cm;                delta         = 2.0 * b2bGap; 
  brickYZ     = (brickY * brickZ);     brickVolume   = (brickYZ * brickX); 
  // Lead Array full length
  PbArrayX    = 360.0*cm;              PbArrayY      = 360.0*cm;          PbArrayZ     = 390.0*cm;
  ArrayXBy2   = 0.5 * PbArrayX;        ArrayYBy2     = 0.5 * PbArrayY;    ArrayZBy2    = 0.5 * PbArrayZ;
  ArraySegX   = PbArrayX / 3.0;        ArraySegY     = PbArrayY / 3.0;    ArraySegZ    = PbArrayZ / 3.0;
  ArraySegX2  = 2.0 * ArraySegX;       ArraySegY2    = 2.0 * ArraySegY;   ArraySegZ2   = 2.0 * ArraySegZ;
  PbArrayXY   = PbArrayX * PbArrayY;   PbArrayVolume = (PbArrayXY * PbArrayZ);
  worldSizeX  = 1.2 * PbArrayX;        worldSizeY    = 1.2 * PbArrayY;    worldSizeZ   = 1.2 * PbArrayZ;

  //----- Here we are using dx, dy and dz as half lengths to minimize computations.
  //---------------------------------------------------------- BoxBottom and BoxTop ------------------------------
  BBDX    =                   PbArrayX - delta;  BBDY =                   PbArrayY - delta;  BBDZ =  ArraySegZ - delta;
  BTDX    =                               BBDX;  BTDY =                               BBDY;  BTDZ =               BBDZ;
  BBCX    =                                0.0;  BBCY =                                0.0;  BBCZ =  -ArrayZBy2 + BBDZ;
  BTCX    =                                0.0;  BTCY =                                0.0;  BTCZ =              -BBCZ;
  BBBTVol =               (BBDX * BBDY * BBDZ);
  //--------------------------------------------------------- RectBottom and RectTop ---------------------------
  RLDX    =                  ArraySegX - delta;  RLDY =                   PbArrayY - delta;  RLDZ =  ArraySegZ - delta;
  RRDX    =                               RLDX;  RRDY =                               RLDY;  RRDZ =               RLDZ;
  RLCX    =            -ArrayXBy2 + 0.5 * RLDX;  RLCY =                                0.0;  RLCZ =                0.0;
  RRCX    =                              -RLCX;  RRCY =                                0.0;  RRCZ =                0.0;
  RLRRVol =               (RLDX * RLDY * RLDZ);
  //---------------------------------------------------------- BoxFront and BoxBack ------------------------------
  BFDX    =                  ArraySegX - delta;  BFDY =            0.5 * ArraySegY - delta;  BFDZ =  ArraySegZ - delta;
  BKDX    =                               BFDX;  BKDY =                               BFDY;  BKDZ =               BFDZ;
  BFCX    =                                0.0;  BFCY =            -ArrayYBy2 + 0.5 * BFDY;  BFCZ =                0.0;
  BKCX    =                                0.0;  BKCY =                              -BFCY;  BKCZ =                0.0;
  BFBKVol =               (BFDX * BFDY * BFDZ);
  //--------------------------------------------------------- BoxCentre -----------------------------------------
  BCDX    =                       BFDX - delta;  BCDY = (PbArrayY - (BFDY + BKDY)) - delta;  BCDZ =       BFDZ - delta;
  BCCX    =                                0.0;  BCCY =                                0.0;  BCCZ =                0.0;
  //--------------------------------------------------------- Hole through BoxFront -----------------------------
  PORTD   = 0.5 * std::min(BFDX, BFDZ) - delta;  PORTL=                       BFDY - delta;
  //---------------------------------------------------------- Cylinder -----------------------------------------
  CYLX    =                                0.0;  CYLY =                                0.0;  CYLZ =                0.0;
  CYIR    =                             0.3*cm;  CYOR =                            10.0*cm;  CYLH =            50.0*cm;
  CBH     =                             0.3*cm;   

  // Compute position Vector for Absorbers, Can and detectors
  posBottomAbsorber    = G4ThreeVector(BBCX, BBCY, BBCZ);
  posTopAbsorber       = G4ThreeVector(BTCX, BTCY, BTCZ);
  posLeftRectAbsorber  = G4ThreeVector(RLCX, RLCY, RLCZ);
  posRightRectAbsorber = G4ThreeVector(RRCX, RRCY, RRCZ);
  posFrontRectAbsorber = G4ThreeVector(BFCX, BFCY, BFCZ);
  posBackRectAbsorber  = G4ThreeVector(BKCX, BKCY, BKCZ);
  posSourceSCan        = G4ThreeVector(BCCX, BCCY, -(BCDZ - (CYLH - 2.0 * CBH)));
}


// Construct
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //Get materials
  Tc                 = G4Material::GetMaterial("G4_Tc");
  Pb                 = G4Material::GetMaterial("G4_Pb");
  steel              = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  Air                = G4Material::GetMaterial("G4_AIR");
  Vac                = G4Material::GetMaterial("G4_Galactic");
  absorberMat        = G4Material::GetMaterial("G4_Pb");
  gapMat             = G4Material::GetMaterial("G4_AIR");
  defaultMat         = G4Material::GetMaterial("G4_AIR");  
 
  
  //------------------------------ declare world
  // World
  auto worldS = new G4Box("World", 0.5 * worldSizeX, 0.5 * worldSizeY, 0.5 * worldSizeZ);  // size
  worldL = new G4LogicalVolume(worldS, defaultMat, "WLogi", 0, 0, 0);
  auto worldP = new G4PVPlacement(0,                 // no rotation
				  G4ThreeVector(),   // at (0,0,0)
				  worldL,           // Logical Volume
				  "worldP",
				  0,                 // Mother volume
				  false,             // no boolean operation
				  0,                 // copy number
				  fCheckOverlaps);   // checking overlaps
   
  // Construction of Lead Absorber and placement of Source Canister
  ConstructAbsorberandCanister();

  // visualization attributes
  G4VisAttributes * Red             = new G4VisAttributes( G4Colour(255/255., 0/255.  , 0/255.   ));
  G4VisAttributes * Yellow          = new G4VisAttributes( G4Colour(255/255., 255/255., 0/255.   ));
  G4VisAttributes * LightBlue       = new G4VisAttributes( G4Colour(0/255.  , 204/255., 204/255. ));
  G4VisAttributes * LightGreen      = new G4VisAttributes( G4Colour(153/255., 255/255., 153/255. ));
  G4VisAttributes * simpleBoxVisAtt = new G4VisAttributes( G4Colour(     1.0,      1.0, 1.0      ));

  worldL->SetVisAttributes (G4VisAttributes::Invisible);
  PbArrayL->SetVisAttributes (G4VisAttributes::Invisible);
  stackL->SetVisAttributes (Yellow);
  simpleBoxVisAtt->SetVisibility(true);
  PbArrayL->SetVisAttributes(simpleBoxVisAtt);

  return worldP;
}

G4VPhysicalVolume* DetectorConstruction::ConstructAbsorberandCanister(){
  // Lead Absorber Array ------ placed at the center which contains all lead bricks 
  auto  PbArrayS      = new G4Box("LeadArray", ArrayXBy2, ArrayYBy2, ArrayZBy2);
  auto  PbArrayL      = new G4LogicalVolume(PbArrayS, defaultMat, "LeadArrayL");     // default Mat (AIR)
  PbArrayP            = new G4PVPlacement(0, 
					  G4ThreeVector(), PbArrayL, "LeadArrayP", 
					  worldL, false, 0, fCheckOverlaps);

  //-------------------- declare BOTTOM and TOP box area ------------------------------->
  auto BigBoxAbsS     = new G4Box("BigBOXS", 0.5 * BBDX, 0.5 * BBDY, 0.5 * BBDZ);
  auto FullAbsorberL  = new G4LogicalVolume(BigBoxAbsS, defaultMat, "FullAbsorberL");
  
  bottomAbsorber      = new G4PVPlacement(0, 
					  posBottomAbsorber, FullAbsorberL, "BottomAbsorber",
					  PbArrayL, false, 0);
  topAbsorber         = new G4PVPlacement(0,
					  posTopAbsorber, fullAbsorberL, "TopAbsorber",
					  PbArrayL, false, 1);
  
  //------------------- delcare LEFT and RIGHT box area ----------------------------->
  auto RectAbsS       = new G4Box("RECTS", 0.5 * RLDX, 0.5 * RLDY, 0.5 * RLDZ);
  auto RectAbsorberL  = new G4LogicalVolume(rectAbsS, defaultMat, "RectAbsL");

  leftAbsorber        = new G4PVPlacement(0,
					  posLeftRectAbsorber, RectAbsorberL, "LeftAbsorber",
					  PbArrayL, false, 0);
  rightAbsorber       = new G4PVPlacement(0,
					  posRightRectAbsorber, RectAbsorberL, "RightAbsorber",
					  PbArrayL, false, 1);

  //------------------- declare FRONT and BACK box area ------------------------------>
  auto SmallBoxAbsS   = new G4Box("SmallBOXS", 0.5 * BFDX, 0.5 * BFDY, 0.5 * BFDZ);
  auto PortS          = new G4Tubs("PORTS", 0.0, 0.5 * PORTD, 0.5 * PORTL, 0.0*deg, 360.0*deg);
  auto RectPortFrontS = new G4UnionSolid("RectPortF", 
					 BoxFS, PortS, RotXP, 
					 G4ThreeVector(0.0, 0.0, 0.0));

  auto RectPortFrontL = new G4LogicalVolume(RectPortfrontS, defaultMat, "RectPortFrontL");
  auto RectBackL      = new G4LogicalVolume(SmallBoxS, defaultMat, "RectBackL");
  frontAbsorber       = new G4PVPlacement(0,
					  posFrontRectAbsorber, RectPortFrontL, "FrontAbsorber",
					  PbArrayL, false, 0);
  backAbsorber        = new G4PVPlacement(0,
					  posBackRectAbsorber, RectBackL, "BackAbsorber",
					  PbArrayL, false, 1);

  //-------------------- declare LLFS contained CYLINDER ----------------------->
  auto CylS           = new G4Tubs("CanS", CYIR, CYOR, 0.5 * CYLH, 0.0*deg, 360.0*deg);
  auto LidS           = new G4Tubs("LidS", 0.0, CYOR, 0.5 * CBH, 0.0*deg, 360.0*deg);
  auto BotS           = new G4Tubs("BotS", 0.0, CYOR, 0.5 * CBH, 0.0*deg, 360.0*deg);
  auto CanLid         = new G4UnionSolid("CanLidU", 
					 CylS, LidS, 0, 
					 G4ThreeVector(0.0, 0.0, 0.5 * (CYLH + CBH)));
  auto LLContS        = new G4UnionSolid("ContBotU", 
					 CanLid, BotS, 0, 
					 G4ThreeVector(0.0, 0.0, -0.5 * (CYLH + CBH)));
  auto LLContL        = new G4LogicalVolume(LLContS, defaultMat, "LLContL");
  cylSrcCan           = new G4PVPlacement(0, 
					  posSourceCan, LLContL, "LLCanP", 
					  PbArrayL, false, 0);  

  return PbArrayP;			  					 
}
					  


// Code for updation of geometries
void DetectorConstruction::UpdateGeometry(){
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  G4RunuManager::getRunManager()->defineWorldVolume(Construct());


}

