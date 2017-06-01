// $Id: DetectorConstruction.cc 2017-05-29 vega $
// 
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
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
  G4double brickLength    = 20.0*cm;
  G4double brickDepth     = 10.0*cm;
  G4double brickHeight    = 5.0*cm;
  G4double b2bGap         = 0.01*cm;
  auto brickXZ            = (brickLength * brickDepth);
  auto brickVolume        = (brickXZ * brickHeight); 
  G4double PbArrayLength  = 150.0*cm;
  G4double PbArrayDepth   = 150.0*cm;
  G4double PbArrayHeight  = 150.0*cm;
  auto ArrayXBy2          = 0.5 * PbArrayLength;
  auto ArrayYBy2          = 0.5 * PbArrayHeight;
  auto ArrayZBy2          = 0.5 * PbArrayDepth;
  auto ArraySegX          = PbArrayLength / 3.0;
  auto ArraySegY          = PbArrayHeight / 3.0;
  auto ArraySegZ          = PbArrayDepth / 3.0;
  auto PbArrayXZ          = PbArrayLength * PbArrayDepth;
  auto PbArrayVolume      = (PbArrayXZ * PbArrayHeight);
  auto nofPblineX         = (int) (PbArrayLength / brickDepth + 0.5);
  auto nofPblineY         = (int) (PbArrayHeight / brickHeight + 0.5);
  auto nofPblineZ         = (int) (PbArrayDepth / brickLength + 0.5);
  auto nofPbXZ            = (int) (PbArrayXZ / brickXZ); 
  auto nofPbVol           = (int) PbArrayVolume / brickVolume; 
  auto nlayers            = (int) nofPbVol / nofPbXZ;
  auto worldSizeX         = 1.2 * PbArrayLength;
  auto worldSizeY         = 1.2 * PbArrayHeight;
  auto worldSizeZ         = 1.2 * PbArrayDepth;

  //----- BoxBottom and BoxTop ------------------------------
  auto BBDX =          PbArrayLength;  auto BBDY =              ArraySegY;  auto BBDZ =                 PbArrayDepth;
  auto BTDX =                   BBDX;  auto BTDY =                   BBDY;  auto BTDZ =                         BBDZ;
  auto BBCX =                    0.0;  auto BBCY = -ArrayYBy2 + ArraySegY;  auto BBCZ =                          0.0;
  auto BTCX =                    0.0;  auto BTCY =  ArrayYBy2 - ArraySegY;  auto BTCZ =                          0.0;
  //------ RectBottom and RectTop ---------------------------
  auto RLDX =              ArraySegX;  auto RLDY =              ArraySegY;  auto RLDZ =                 PbArrayDepth;
  auto RRDX =                   RLDX;  auto RRDY =                   RLDY;  auto RRDZ =                         RLDZ;
  auto RLCX = -ArrayXBy2 + ArraySegX;  auto RLCY =                    0.0;  auto RLCZ =                          0.0;
  auto RRCX =  ArrayXBy2 - ArraySegX;  auto RRCY =                    0.0;  auto RRCZ =                          0.0;
  //----- BoxFront and BoxBack ------------------------------
  auto BFDX =              ArraySegX;  auto BFDY =              ArraySegY;  auto BFDZ =              0.5 * ArraySegZ;
  auto BKDX =                   BFDX;  auto BKDY =                   BFDY;  auto BKDZ =                         BFDZ;
  auto BFCX =                    0.0;  auto BFCY =                    0.0;  auto BFCZ =            -ArrayZBy2 + BFDZ;
  auto BKCX =                    0.0;  auto BKCY =                    0.0;  auto BKCZ =             ArrayZBy2 - BFDZ;
  //----- Cylinder -----------------------------------------
  auto CYLX =                    0.0;  auto CYLY =                    0.0;  auto CYLZ =                          0.0;
  auto CYIR =                 0.3*cm;  auto CYOR =                10.0*cm;  auto CYLH =                      40.0*cm;
  auto CBH  =                 0.3*cm;   

  //Get materials
  auto Pb                 = G4Material::GetMaterial("G4_Pb");
  auto steel              = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto Air                = G4Material::GetMaterial("G4_AIR");
  auto Vac                = G4Material::GetMaterial("G4_Galactic");
  auto absorberMat        = G4Material::GetMaterial("G4_Pb");
  auto gapMat             = G4Material::GetMaterial("G4_AIR");
  auto defaultMat         = G4Material::GetMaterial("G4_AIR");  

   
  G4int i = 0, j = 0, k = 0, tmp1 = 0, tmp2 = 0, tmp3 = 0;
  G4double xnew= 0., ynew = 0., znew = 0.;
  G4double x0= -0.5 * PbArrayLength, y0 = -0.5 * PbArrayHeight, z0 = -0.5*PbArrayDepth;

  //------------------------------ declare world
  // World
  auto worldS = new G4Box("World", 0.5*worldSizeX, 0.5*worldSizeY, 0.5*worldSizeZ);  // size
  auto worldLV = new G4LogicalVolume(worldS, defaultMat, "WLogi", 0, 0, 0);
  auto worldPV = new G4PVPlacement(0,                 // no rotation
				   G4ThreeVector(),   // at (0,0,0)
				   worldLV,           // Logical Volume
				   "wPhys",
				   0,                 // Mother volume
				   false,             // no boolean operation
				   0,                 // copy number
				   fCheckOverlaps);   // checking overlaps
   
  //-------------------- declare Lead Array
  // Lead Absorber Array ------ placed at the center which contains all lead bricks 
  auto  PbArrayS = new G4Box("LeadArray", 0.5 * PbArrayLength, 0.5 * PbArrayHeight, 0.5 * PbArrayDepth);
  auto  PbArrayLV = new G4LogicalVolume(PbArrayS, defaultMat, "LeadArrayL");             // default Material given (AIR)
  new G4PVPlacement(0, G4ThreeVector(), PbArrayLV, "LeadArrayP", worldLV, false, 0, fCheckOverlaps);

  //-------------------- declare BOTTOM box area ------------------------------->
  auto BoxBS = new G4Box("BOXBS", BBDX, BBDY, BBDZ);
  auto BoxBL = new G4LogicalVolume(BoxBS, defaultMat, "BoxBackL");
  new G4PVPlacement(0, G4Threevector(BBCX, BBCY, BBCZ), "BoxBP", PbArrayLV, false, 0, fCheckOverlaps);
  //-------------------- declare TOP box area ---------------------------------->
  auto BoxTS = new G4Box("BOXTS", 0.5 *BTDX, 0.5 * BTDY, 0.5 * BTDZ);
  auto BoxTL = new G4LogicalVolume(BoxTS, defaultMat, "BoxTopL");
  new  G4PVPlacement(0, G4ThreeVector(BTCX, BTCY, BTCZ), "BoxTP", PbArrayLV, false, 0, fCheckOverlaps);
  //-------------------- declare LEFT rectangular paralleopipped area ---------->
  auto RectLS = new G4Box("RECTLS", 0.5 * RLDX, 0.5 * RLDY, 0.5 * RLDZ);
  auto RectLL = new G4LogicalVolume(RectLS, defaultMat, "RectLeftL");
  new  G4PVPlacement(0, G4ThreeVector(RLCX, RLCY, RLCZ), "RectLP", PbArrayLV, false, 0, fCheckOverlaps);
  //-------------------- delcare RIGHT rectangular parallelopipped area -------->
  auto RectRS = new G4Box("RECTRS", 0.5 * RRDX, 0.5 * RRDY, 0.5 * RRDZ);
  auto RectRL = new G4LogicalVolume(RectRS, defaultMat, "RectRightL");
  new G4PVPlacement(0, G4ThreeVector(RRCX, RRCY, RRCZ), "RectRP", PbArrayLV, false, 0, fCheckOverlaps);
  //-------------------- declare FRONT box area -------------------------------->
  auto BoxFS = new G4Box("BOXFS", 0.5 * BFDX, 0.5 * BFDY, 0.5 * BFDZ);
  auto BoxFL = new G4LogicalVolume(BoxFS, defaultMat, "BoxFrontL");
  new G4PVPlacement(0, G4ThreeVector(BFCX, BFCY, BFCZ), "BoxFP", PbArrayLV, false, 0, fCheckOverlaps);
  //-------------------- declare BACK box area --------------------------------->
  auto BoxKS = new G4Box("BOXKS", 0.5 * BKDX, 0.5 * BKDY, 0.5 * BKDZ);
  auto BoxKL = new G4LogicalVolume(BoxKS, defaultMat, "BoxBackL");
  new G4PVPlacement(0, G4ThreeVector(BKCX, BKCY, BKCZ), "BoxKP", PbArrayLV, false, 0, fCheckOverlaps);
  //-------------------- declare LLFS contained CYLINDER ----------------------->
  auto CylS = new G4Tubs("CanS", CYIR, CYOR, 0.5 * (CYLH - 2.0 * CBH), 0.0*deg, 360.0*deg);
  auto LidS = new G4Tubs("LidS", 0.0, CYOR, 0.5 * CBH, 0.0*deg, 360.0*deg);
  auto BotS = new G4Tubs("BotS", 0.0, CYOR, 0.5 * CBH, 0.0*deg, 360.0*deg);
  //---- Move Lid to top of the Can and make Union Solid
  auto CanLid = new G4UnionSolid(&CylS, &LidS, 0, G4ThreeVector(0.0, 0.0, 0.5 * (CYLH - 2.0 * CBH)));
  auto LLContS = new G4UnionSolid(&CanLis, &BotS, 0, G4ThreeVector(0.0, 0.0, -0.5 * (CYLH - 2.0 * CBH)));
  


  //--------------------- declare stacks where stacks of bricks are laid
  // Stack size brick size + 1/2 of gap between two consecutive bricks  << Arrays of stacks >>
  auto stackS = new G4Box("Stack", 0.5*brickDepth + 0.5*b2bGap,
			  0.5*brickHeight + 0.5*b2bGap,
			  0.5*brickLength + 0.5*b2bGap);
  auto stackL = new G4LogicalVolume(stackS, defaultMat, "stackL");
  i = j = k = 1;
  for (int Cnt = 0; Cnt < nofPbVol; Cnt++){
    G4double delX = i * (brickDepth + 0.5 * b2bGap);
    G4double delY = (j - 1) * (brickHeight + 0.5 * b2bGap);
    G4double delZ = (k - 1) * (brickLength + 0.5 * b2bGap);
    xnew = x0 + delX;    ynew = y0 + delY;   znew = z0 + delZ;
    new G4PVPlacement(0, G4ThreeVector(xnew, ynew, znew), stackL, "stack", PbArrayLV, false, Cnt, fCheckOverlaps);
    i = ((Cnt % nofPblineX) == 0) ? 1 : i;
    j = ((Cnt % nofPblineX) == 0) ? ++j : j;
    k = ((Cnt % nofPbXZ) == 0) ? ++k : k;
  }
	 
  //-------------------- declare lead bricks
  //Lead bricks ------  placed in array within the Lead Array object  << single objects >>
  auto brickS = new G4Box("brickS", 0.5 * brickDepth, 0.5 * brickHeight, 0.5 * brickLength);
  auto brickLV = new G4LogicalVolume(brickS, Pb, "brickL");                                              // brick of Lead (Pb)
  auto brickP = new G4PVPlacement(0, G4ThreeVector(xnew, ynew, znew), brickLV,"brick", stackL, false,i,fCheckOverlaps);  


  // visualization attributes

  G4VisAttributes * Red             = new G4VisAttributes( G4Colour(255/255., 0/255.  , 0/255.   ));
  G4VisAttributes * Yellow          = new G4VisAttributes( G4Colour(255/255., 255/255., 0/255.   ));
  G4VisAttributes * LightBlue       = new G4VisAttributes( G4Colour(0/255.  , 204/255., 204/255. ));
  G4VisAttributes * LightGreen      = new G4VisAttributes( G4Colour(153/255., 255/255., 153/255. ));
  G4VisAttributes * simpleBoxVisAtt = new G4VisAttributes( G4Colour(     1.0,      1.0, 1.0      ));

  worldLV->SetVisAttributes (G4VisAttributes::Invisible);
  PbArrayLV->SetVisAttributes (G4VisAttributes::Invisible);
  stackL->SetVisAttributes (Yellow);

  
  simpleBoxVisAtt->SetVisibility(true);
  PbArrayLV->SetVisAttributes(simpleBoxVisAtt);

  return worldPV;
}




