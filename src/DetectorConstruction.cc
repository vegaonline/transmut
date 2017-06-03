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


//
G4ThreadLocal

// The Constructor
DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction(),
					       fAbsorberPV(0),
					       fTargetPV(0),
					       fCheckOverlaps(true)
{
}

// The Destructor
DetectorConstruction::~DetectorConstruction()
{
}

// Construct
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define Materials
  DefineMaterials();

  // Define Volumes
  return DefineVolumes();
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

// Define Volumes
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry Parameters
  // brick full length
  G4double brickX         = 20.0*cm;
  G4double brickY         = 10.0*cm;
  G4double brickZ         = 5.0*cm;
  auto HalfBrickX         = 0.5 * brickX;
  auto HalfBrickY         = 0.5 * brickY;
  auto HalfBrickZ         = 0.5 * brickZ;
  G4double b2bGap         = 0.01*cm;
  auto brickYZ            = (brickY * brickZ);
  auto brickVolume        = (brickYZ * brickX); 
  // Lead Array full length
  G4double PbArrayX       = 360.0*cm;
  G4double PbArrayY       = 360.0*cm;
  G4double PbArrayZ       = 390.0*cm;

  auto ArrayXBy2          = 0.5 * PbArrayX;
  auto ArrayYBy2          = 0.5 * PbArrayY;
  auto ArrayZBy2          = 0.5 * PbArrayZ;
  auto ArraySegX          = PbArrayX / 3.0;
  auto ArraySegY          = PbArrayY / 3.0;
  auto ArraySegZ          = PbArrayZ / 3.0;
  auto ArraySegX2         = 2.0 * ArraySegX;
  auto ArraySegY2         = 2.0 * ArraySegY;
  auto ArraySegZ2         = 2.0 * ArraySegZ;
  auto PbArrayXY          = PbArrayX * PbArrayY;
  auto PbArrayVolume      = (PbArrayXY * PbArrayZ);
  auto worldSizeX         = 1.2 * PbArrayX;
  auto worldSizeY         = 1.2 * PbArrayY;
  auto worldSizeZ         = 1.2 * PbArrayZ;

  //----- Here we are using dx, dy and dz as half lengths to minimize computations.
  //-------------------------------------------- BoxBottom and BoxTop ------------------------------
  auto BBDX    =                   PbArrayX;  auto BBDY =                    PbArrayY;  auto BBDZ =          ArraySegZ;
  auto BTDX    =                       BBDX;  auto BTDY =                        BBDY;  auto BTDZ =               BBDZ;
  auto BBCX    =                        0.0;  auto BBCY =                         0.0;  auto BBCZ =  -ArrayZBy2 + BBDZ;
  auto BTCX    =                        0.0;  auto BTCY =                         0.0;  auto BTCZ =              -BBCZ;
  auto BBBTVol =       (BBDX * BBDY * BBDZ);
  //-------------------------------------------- RectBottom and RectTop ---------------------------
  auto RLDX    =                  ArraySegX;  auto RLDY =                    PbArrayY;  auto RLDZ =          ArraySegZ;
  auto RRDX    =                       RLDX;  auto RRDY =                        RLDY;  auto RRDZ =               RLDZ;
  auto RLCX    =    -ArrayXBy2 + 0.5 * RLDX;  auto RLCY =                         0.0;  auto RLCZ =                0.0;
  auto RRCX    =                      -RLCX;  auto RRCY =                         0.0;  auto RRCZ =                0.0;
  auto RLRRVol =       (RLDX * RLDY * RLDZ);
  //-------------------------------------------- BoxFront and BoxBack ------------------------------
  auto BFDX    =                  ArraySegX;  auto BFDY =             0.5 * ArraySegY;  auto BFDZ =          ArraySegZ;
  auto BKDX    =                       BFDX;  auto BKDY =                        BFDY;  auto BKDZ =               BFDZ;
  auto BFCX    =                        0.0;  auto BFCY =     -ArrayYBy2 + 0.5 * BFDY;  auto BFCZ =                0.0;
  auto BKCX    =                        0.0;  auto BKCY =                       -BFCY;  auto BKCZ =                0.0;
  auto BFBKVol =       (BFDX * BFDY * BFDZ);
  //-------------------------------------------- BoxCentre -----------------------------------------
  auto BCDX    =                       BFDX;  auto BCDY =  (PbArrayY - (BFDY + BKDY));  auto BCDZ =               BFDZ;
  auto BCCX    =                        0.0;  auto BCCY =                         0.0;  auto BCCZ =                0.0;
  //-------------------------------------------- Hole through BoxFront -----------------------------
  auto PORTD   = 0.5 * std::min(BFDX, BFDZ);  auto PORTL=                        BFDY;
  //-------------------------------------------- Cylinder -----------------------------------------
  auto CYLX    =                        0.0;  auto CYLY =                         0.0;  auto CYLZ =                0.0;
  auto CYIR    =                     0.3*cm;  auto CYOR =                     10.0*cm;  auto CYLH =            50.0*cm;
  auto CBH     =                     0.3*cm;   

  //Get materials
  auto Tc                 = G4Material::GetMaterial("G4_Tc");
  auto Pb                 = G4Material::GetMaterial("G4_Pb");
  auto steel              = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto Air                = G4Material::GetMaterial("G4_AIR");
  auto Vac                = G4Material::GetMaterial("G4_Galactic");
  auto absorberMat        = G4Material::GetMaterial("G4_Pb");
  auto gapMat             = G4Material::GetMaterial("G4_AIR");
  auto defaultMat         = G4Material::GetMaterial("G4_AIR");  

   
  G4RotationMatrix *RotXP = new G4RotationMatrix();  
  G4RotationMatrix *RotYP = new G4RotationMatrix();
  G4RotationMatrix *RotZP = new G4RotationMatrix();
  G4RotationMatrix *RotXM = new G4RotationMatrix();
  G4RotationMatrix *RotYM = new G4RotationMatrix();
  G4RotationMatrix *RotZM = new G4RotationMatrix();
  RotXP->rotateX( 90*deg);   RotYP->rotateY( 90*deg);  RotZP->rotateZ( 90*deg);
  RotXM->rotateX(-90*deg);   RotYM->rotateY(-90*deg);  RotZM->rotateZ(-90*deg);

  G4int i = 0;
  G4int j = 0;
  G4int k = 0;
  G4double xnew= 0.0;
  G4double ynew = 0.0;
  G4double znew = 0.0;
  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  G4int nPbX = 0;
  G4int nPbY = 0;
  G4int nPbZ = 0;
  G4int nPb = 0;
 

  //------------------------------ declare world
  // World
  auto worldS = new G4Box("World", 0.5 * worldSizeX, 0.5 * worldSizeY, 0.5 * worldSizeZ);  // size
  auto worldL = new G4LogicalVolume(worldS, defaultMat, "WLogi", 0, 0, 0);
  auto worldP = new G4PVPlacement(0,                 // no rotation
				   G4ThreeVector(),   // at (0,0,0)
				   worldL,           // Logical Volume
				   "worldP",
				   0,                 // Mother volume
				   false,             // no boolean operation
				   0,                 // copy number
				   fCheckOverlaps);   // checking overlaps
   
  //-------------------- declare Lead Array
  // Lead Absorber Array ------ placed at the center which contains all lead bricks 
  auto  PbArrayS = new G4Box("LeadArray", ArrayXBy2, ArrayYBy2, ArrayZBy2);
  auto  PbArrayL = new G4LogicalVolume(PbArrayS, defaultMat, "LeadArrayL");             // default Material given (AIR)
  new G4PVPlacement(0, G4ThreeVector(), PbArrayL, "LeadArrayP", worldL, false, 0, fCheckOverlaps);

  //-------------------- declare BOTTOM box area ------------------------------->
  auto BoxBS = new G4Box("BOXBS", 0.5 * BBDX, 0.5 * BBDY, 0.5 * BBDZ);
  auto BoxBL = new G4LogicalVolume(BoxBS, defaultMat, "BoxBackL");
  new G4PVPlacement(0, G4ThreeVector(BBCX, BBCY, BBCZ), BoxBL, "BoxBP", PbArrayL, false, 0, fCheckOverlaps);

  //-------------------- declare LEFT rectangular paralleopipped area ---------->
  auto RectLS = new G4Box("RECTLS", 0.5 * RLDX, 0.5 * RLDY, 0.5 * RLDZ);
  auto RectLL = new G4LogicalVolume(RectLS, defaultMat, "RectLeftL");
  new  G4PVPlacement(0, G4ThreeVector(RLCX, RLCY, RLCZ), RectLL, "RectLP", PbArrayL, false, 0, fCheckOverlaps);

  //-------------------- declare FRONT box with a hole area -------------------------------->
  auto BoxFS = new G4Box("BOXFS", 0.5 * BFDX, 0.5 * BFDY, 0.5 * BFDZ);
  auto BoxFL = new G4LogicalVolume(BoxFS, defaultMat, "BoxFL");

//-------- Now place hole inside stacked BOXF for creating port ------------->
  auto PortS = new G4Tubs("PORTS", 0.0, 0.5 * PORTD, 0.5 * PORTL, 0.0*deg, 360.0*deg);
  auto BFPortS = new G4UnionSolid("BFPortU", BoxFS, PortS, RotXP, G4ThreeVector(0.0, 0.0, 0.0));
  auto BFPortL = new G4LogicalVolume(BFPortS, defaultMat, "BFPortL");
  new G4PVPlacement(0, G4ThreeVector(BFCX, BFCY, BFCZ), BFPortL, "BFPortP", PbArrayL, false, 0, fCheckOverlaps);

  //-------------------- declare BACK box area --------------------------------->
  auto BoxKS = new G4Box("BOXKS", 0.5 * BKDX, 0.5 * BKDY, 0.5 * BKDZ);
  auto BoxKL = new G4LogicalVolume(BoxKS, defaultMat, "BoxBackL");
  new G4PVPlacement(0, G4ThreeVector(BKCX, BKCY, BKCZ), BoxKL, "BoxKP", PbArrayL, false, 0, fCheckOverlaps);

  //-------------------- declare Central box for placing LLFS and detectors ---->
  auto BoxCS = new G4Box("BOXCS", 0.5 * BCDX, 0.5 * BCDY, 0.5 * BCDZ);
  auto BoxCL = new G4LogicalVolume(BoxCS, defaultMat, "BoxCL");
  new G4PVPlacement(0, G4ThreeVector(BCCX, BCCY, BCCZ), BoxCL, "BoxCP", PbArrayL, false, 0, fCheckOverlaps);

  //-------------------- delcare RIGHT rectangular parallelopipped area -------->
  auto RectRS = new G4Box("RECTRS", 0.5 * RRDX, 0.5 * RRDY, 0.5 * RRDZ);
  auto RectRL = new G4LogicalVolume(RectRS, defaultMat, "RectRightL");
  new G4PVPlacement(0, G4ThreeVector(RRCX, RRCY, RRCZ), RectRL, "RectRP", PbArrayL, false, 0, fCheckOverlaps);

  //-------------------- declare TOP box area ---------------------------------->
  auto BoxTS = new G4Box("BOXTS", 0.5 * BTDX, 0.5 * BTDY, 0.5 * BTDZ);
  auto BoxTL = new G4LogicalVolume(BoxTS, defaultMat, "BoxTopL");
  new  G4PVPlacement(0, G4ThreeVector(BTCX, BTCY, BTCZ), BoxTL, "BoxTP", PbArrayL, false, 0, fCheckOverlaps);

  //-------------------- declare LLFS contained CYLINDER ----------------------->
  auto CylS = new G4Tubs("CanS", CYIR, CYOR, 0.5 * (CYLH - 2.0 * CBH), 0.0*deg, 360.0*deg);
  auto LidS = new G4Tubs("LidS", 0.0, CYOR, 0.5 * CBH, 0.0*deg, 360.0*deg);
  auto BotS = new G4Tubs("BotS", 0.0, CYOR, 0.5 * CBH, 0.0*deg, 360.0*deg);
  auto CanLid = new G4UnionSolid("CanLidU", CylS, LidS, 0, G4ThreeVector(0.0, 0.0, 0.5 * (CYLH - 2.0 * CBH)));
  auto LLContS = new G4UnionSolid("ContBotU", CanLid, BotS, 0, G4ThreeVector(0.0, 0.0, -0.5 * (CYLH - 2.0 * CBH)));
  auto LLContL = new G4LogicalVolume(LLContS, defaultMat, "LLContL");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -4.7*cm), LLContL, "LLCanP", BoxCL, false, 0, fCheckOverlaps);  
	 
  //-------------------- declare lead bricks
  //Lead bricks ------  placed in array within the Lead Array object  << single objects >>
  auto brickS = new G4Box("brickS", HalfBrickX, HalfBrickY, HalfBrickZ);
  auto brickL = new G4LogicalVolume(brickS, Pb, "brickL");                                 // brick of Lead (Pb)

  //--------------------- declare stacks where stacks of bricks are laid
  // Stack size brick size + 1/2 of gap between two consecutive bricks  << Arrays of stacks >>
  auto stackS = new G4Box("StackS", HalfBrickX, HalfBrickY, HalfBrickZ);
  auto stackL = new G4LogicalVolume(stackS, defaultMat, "stackL");
  
  //------ STACKS :: BB + BT ------>
  i = j = k = 1;
  x0 = -ArrayXBy2, y0 = -ArrayYBy2, z0 = -ArrayZBy2;
  nPbX = (int)(BBDX / brickX);             // brickX along BBX
  nPbY = (int)(BBDY / brickY);             // brickY along BBY
  nPbZ = (int)(BBDZ / brickZ);             // brickZ along BBZ
  nPb  = (nPbX * nPbY * nPbZ);       // Total bricks required
  for (int Cnt = 0; Cnt < nPb; Cnt++){
    xnew = x0 +       i * brickX;
    ynew = y0 + (j - 1) * brickY;
    znew = z0 + (k - 1) * brickZ;    
    new G4PVPlacement(0, G4ThreeVector(xnew, ynew, znew), stackL, "stackBBP", BoxBL, false, 0, fCheckOverlaps);
    new G4PVPlacement(0, G4ThreeVector(xnew, ynew, znew + ArraySegZ2), stackL, "stackBTP", BoxTL, false, 0, fCheckOverlaps);
    x0 = xnew;  y0 = ynew;  z0 = znew;
    i = ((Cnt % nPbX) == 0) ?   1 : i;
    j = ((Cnt % nPbX) == 0) ? ++j : j;
    k = ((Cnt % nPbZ) == 0) ? ++k : k;
  }
  //------ STACKS :: RL + RT ------>
  i = j = k = 1;
  x0 = -ArrayXBy2, y0 = -ArrayYBy2, z0 = -0.5 * RLDZ;
  nPbX = (int)(RLDX / brickY);             // brickY along RLX RRX
  nPbY = (int)(RLDY / brickX);             // brickY along RLY RRY
  nPbZ = (int)(RLDZ / brickZ);             // brickZ along RLZ RRZ
  nPb  = (nPbX * nPbY * nPbZ);       // Total bricks required
  for (int Cnt = 0; Cnt < nPb; Cnt++){
    xnew = x0 +       i * brickY;
    ynew = y0 + (j - 1) * brickX;
    znew = z0 + (k - 1) * brickZ;    
    new G4PVPlacement(RotZP, G4ThreeVector(xnew, ynew, znew)             , stackL, "stackRLP", RectLL, false, 0, fCheckOverlaps);
    new G4PVPlacement(RotZP, G4ThreeVector(xnew + ArraySegX2, ynew, znew), stackL, "stackRRP", RectRL, false, 0, fCheckOverlaps);
    x0 = xnew;  y0 = ynew;  z0 = znew;
    i = ((Cnt % nPbX) == 0) ?   1 : i;
    j = ((Cnt % nPbX) == 0) ? ++j : j;
    k = ((Cnt % nPbZ) == 0) ? ++k : k;
  }
  //-------- STACKS :: BFPortL + BK ------>
  i = j = k = 1;
  x0 = -0.5 * BFDX, y0 = -0.5 * BFDY, z0 = -0.5 * BFDZ;
  nPbX = (int)(BFDX / brickX);             // brickX along BFX BKX
  nPbY = (int)(BFDY / brickY);             // brickY along BFY BKY
  nPbZ = (int)(BFDZ / brickZ);             // brickZ along BFZ BKZ
  nPb  = (int)( nPbX * nPbY * nPbZ);       // Total bricks required
  for (int Cnt = 0; Cnt < nPb; Cnt++){
    xnew = x0 +       i * brickY;
    ynew = y0 + (j - 1) * brickX;
    znew = z0 + (k - 1) * brickZ;    
    new G4PVPlacement(RotZP, G4ThreeVector(xnew, ynew, znew)             , stackL, "stackBFP", BFPortL, false, 0, fCheckOverlaps);
    new G4PVPlacement(RotZP, G4ThreeVector(xnew, ynew + ArraySegY2, znew), stackL, "stackBKP", BoxKL, false, 0, fCheckOverlaps);
    x0 = xnew;  y0 = ynew;  z0 = znew;
    i = ((Cnt % nPbX) == 0) ?   1 : i;
    j = ((Cnt % nPbX) == 0) ? ++j : j;
    k = ((Cnt % nPbZ) == 0) ? ++k : k;
  }  
  


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




