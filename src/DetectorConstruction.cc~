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
    G4Nistmanager* nistManager = G4NistManager::Instance();

    // Select Lead 100% pure
    G4Material* Pb             = nistManager->FindOrBuildMaterial("G4_Pb);
    G4Material* Be             = nistManager->FindOrBuildMaterial("G4_Be");
    G4Material* polyFilm       = nistManager->FindOrBuidlMaterial("G4_KAPTON");
    G4Material* steel          = nistManager->FindOrBuildMaterial("G4_Stainless-Steel");
    G4Material* Air            = nistManager->FindOrBuildMaterial("G4_AIR");
    G4Material* Argongas       = nistManager->FindOrBuildMaterial("G4_Ar");
    G4Material* Vac            = nistmanager->FindOrBuildMaterial("G4_Galactic");
    G4Material* defaultMat     = nistManager->FindOrBuildMaterial("G4_Galactic");        

    // Print the material synthesized
    G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
    G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

// Defien Volumes
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
    // Geometry Parameters
    G4double brickLength    = 20.0*cm;
    G4double brickDepth     = 10.0*cm;
    G4double brickHeight    = 5.0*cm;
    G4double brickXY        = (brickLength * brickheight);
    G4double brickVolume    = (brickXY * brickDepth); 
    G4double PbArrayLength  = 150.0*cm;
    G4double PbArrayDepth   = 150.0*cm;
    G4double PbArrayHeight  = 150.0*cm;
    G4double PbArrayXY      = PbArrayLength * PbArrayHeight;
    G4double PbArrayVolume  = (PbArrayXY * PbArrayDepth);
    G4int nofPblineX        = (int) (PbArrayLength / brickLength + 0.5);
    G4int nofPblineY        = (int) (PbArrayHeight / brickHeight + 0.5);
    G4int nofPbXY           = (PbArrayXY / brickXY); 
    G4int nofPbVol          = PbArrayVolume / brickVolume; 
    G4int nofPblineZ        = nofPbVol / nofPbXY;
    G4int nlayers           = nofPbVol / nofPblineX;
    G4double worldSizeXY    = 1.2 * PbArrayXY;
    G4double worldSizeZ     = 1.2 * PbArrayDepth;

   
    unsigned i = 0, j = 0, k = 0, tmp1 = 0, tmp2 = 0;
    G4double xnew= 0., ynew = 0., znew = 0.;

    // World
    G4VSolid* worldS = new G4Box("World", 0.5*worldSizeXY, 0.5*worldSizeXY, 0.5*worldSizeXY);  // size
    G4VLogicalVolume* worldLV = new G4LogicalVolume(worldS, defaultMat, "WLog", 0, 0, 0);
    G4VPhysicalVolume* worldPV = new G4PVPlacement(0,                 // no rotation
                                                   G4ThreeVector(),   // at (0,0,0)
                                                   worldLV,           // Logical Volume
                                                   "wPhys",
                                                   0,                 // Mother volume
                                                   false,             // no boolean operation
                                                   0,                 // copy number
                                                   fCheckOverlaps);   // checking overlaps
   
    // Lead Absorber
    G4VSolid* PbArray = new G4Box("LeadArray", 0.5 * PbArrayXy, 0.5 * PbArrayXY, 0.5 * PbArraydepth);
    G4VLogicalVolume* PbArrayLV = new G4LogicalVolume(PbArray, Pb, "LeadArray");
    for (i = 0; i < nofPbVol; i++) 
    {
        j = (j < nofPblineX) ? J++ : 0;
        tmp1 = j;       tmp2 = i / nofPblineX;      tmp3 = i / nofPbXY;
        xnew = tmp1 * brickLength;     ynew = tmp2 * brickHeight;    znew = tmp3 * brickDepth;
        new G4PVPlacement(0, G4ThreeVectors(xnew, ynew, znew), PbArrayLV, 
                          "PbArrayPhysical", worldLV, false,i,fCheckOverlaps);
    }


    // visualization attributes
    worldLV->SetVisAttributes (G4VisAttributes::Invisible);
    G4VisAttributes* simpleBoxVisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
    simpleBoxVisAtt->SetVisibility(true);
    PbArrayLV->SetVisAttributes(simpleBoxVisAtt);

    return worldPV;
}




