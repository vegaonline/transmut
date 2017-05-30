// $Id: DetectorConstruction.cc 2014-11-26 17:08:44Z ABHIJIT $
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
  auto brickXZ            = (brickLength * brickDepth);
  auto brickVolume        = (brickXZ * brickHeight); 
  G4double PbArrayLength  = 150.0*cm;
  G4double PbArrayDepth   = 150.0*cm;
  G4double PbArrayHeight  = 150.0*cm;
  auto PbArrayXZ          = PbArrayLength * PbArrayDepth;
  auto PbArrayVolume      = (PbArrayXZ * PbArrayHeight);
  auto nofPblineX         = (int) (PbArrayLength / brickDepth + 0.5);
  auto nofPblineY         = (int) (PbArrayHeight / brickHeight + 0.5);
  auto nofPblineZ         = (int) (PbArrayDepth / brickLength + 0.5);
  auto nofPbXZ            = (PbArrayXZ / brickXZ); 
  auto nofPbVol           = PbArrayVolume / brickVolume; 
  auto nlayers            = nofPbVol / nofPbXZ;
  auto worldSizeX         = 1.2 * PbArrayLength;
  auto worldSizeY         = 1.2 * PbArrayHeight;
  auto worldSizeZ         = 1.2 * PbArrayDepth;

  //Get materials
  auto Pb                 = G4Material::GetMaterial("G4_Pb");
  auto steel              = G4Material::GetMaterial("G4_STAINLESS-STEEL");
  auto Air                = G4Material::GetMaterial("G4_AIR");
  auto Vac                = G4Material::GetMaterial("G4_Galactic");
  auto absorberMat        = G4Material::GetMaterial("G4_Pb");
  auto gapMat             = G4Material::GetMaterial("Air");
  auto defaultMat         = G4Material::GetMaterial("Air");  

   
  unsigned i = 0, j = 0, k = 0, tmp1 = 0, tmp2 = 0, tmp3 = 0;
  G4double xnew= 0., ynew = 0., znew = 0.;
  G4double x0= -0.5 * PbArrayLength, y0 = -0.5 * PbArrayDepth, z0 = -0.5*PbArrayHeight;

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
  // Lead Absorber Array ------ placed at the center which contains all lead bricks << SINGLE OBJECT >>
  auto  PbArrayS = new G4Box("LeadArrayS", 0.5 * PbArrayLength, 0.5 * PbArrayHeight, 0.5 * PbArrayDepth);
  auto  PbArrayLV = new G4LogicalVolume(PbArrayS, defaultMat, "LeadArrayL");             // default Material given (AIR)
  new G4PVPlacement(0, G4ThreeVector(), PbArrayLV, "LeadArray", worldLV, false, 0, fCheckOverlaps);


  //-------------------- declare lead bricks
  //Lead bricks ------  placed in array within the Lead Array object  << Arrays of objects >>
  auto brickS = new G4Box("brickS", 0.5 * brickDepth, 0.5 * brickHeight, 0.5 * brickLength);
  auto brickLV = new G4LogicalVolume(brickS, Pb, "brickL");                                              // brick of Lead (Pb)
  auto brickP = new G4PVPlacement(0, G4ThreeVector(xnew, ynew, znew), PbArrayLV, 
				  "PbArrayPhys", worldLV, false,i,fCheckOverlaps);  

	
  for (i = 0; i < nofPbVol; i++) {
    j = (j < nofPblineX) ? j++ : 0;
    tmp1 = j;       tmp2 = i / nofPblineX;      tmp3 = i / nofPbXZ;
    xnew = x0 + tmp1 * brickDepth;     ynew = y0 + tmp2 * brickLength;    znew = z0 + tmp3 * brickHeight;
    auto brickP = new G4PVPlacement(0, G4ThreeVector(xnew, ynew, znew), PbArrayLV, 
				    "PbArrayPhys", worldLV, false,i,fCheckOverlaps);
  }


  // visualization attributes

  G4VisAttributes * Red        = new G4VisAttributes( G4Colour(255/255. ,0/255.   ,0/255.   ));
  G4VisAttributes * Yellow     = new G4VisAttributes( G4Colour(255/255. ,255/255. ,0/255.   ));
  G4VisAttributes * LightBlue  = new G4VisAttributes( G4Colour(0/255.   ,204/255. ,204/255. ));
  G4VisAttributes * LightGreen = new G4VisAttributes( G4Colour(153/255. ,255/255. ,153/255. ));
  G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));

  worldLV->SetVisAttributes (G4VisAttributes::Invisible);
  PbArrayLV->SetVisAttributes (G4VisAttributes::Invisible);
  brickLV->SetVisAttributes (Yellow);

  
  simpleBoxVisAtt->SetVisibility(true);
  PbArrayLV->SetVisAttributes(simpleBoxVisAtt);

  return worldPV;
}




