// $Id: DetectorConstruction.hh 2017-05-29 vega $
// // 
// /// \file DetectorConstruction.hh
// /// \brief Definition of the DetectorConstruction class
//
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class G4FieldManager;
class G4GlobalMagFieldMessenger;

// The detector is a box of High purity Lead 99.999%
// Parameters to control geometry of the calorimetry
// Lead chamber dimension
// Beam Parameters

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

public:
  G4VPhysicalVolume* Construct();       // construct geometry of the setup
  //        virtual void ConstructSDandField();
  void UpdateGeometry();


  const G4VPhysicalVolume* GetAbsorberPV() const;
  const G4VPhysicalVolume* GetTargetPV() const;

private:
  // Methods
  //
  void DefineMaterials();           // Define materials to be filled in geometries
  void ComputeParams();             // Define dimensions and position coordinates for geometries
  G4VPhysicalVolume* ConstructAbsorberandCanister();    // make geometries

  // Data members
  //
  // static G4ThreadLocal G4GlobalMagFieldMessenger* fMagFieldMessenger;
  G4VPhysicalVolume* fAbsorberPV;           // the absorber physical volume
  G4VPhysicalVolume* fTargetPV;             // the LLFS target to produce neutron

  G4bool fCheckOverlaps;                   // check volumes overlap

private:

  //dummies
  G4int i = 0;
  G4int j = 0;
  G4int k = 0;
  G4double xnew= 0.0;
  G4double ynew = 0.0;
  G4double znew = 0.0;
  G4double x0 = 0.0;
  G4double y0 = 0.0;
  G4double z0 = 0.0;
  int nPbX = 0;
  int nPbY = 0;
  int nPbZ = 0;
  int nPb = 0;

  //materials
  G4Material* Tc;
  G4Material*  Pb;
  G4Material* steel;
  G4Material* Air;
  G4Material* Vac;
  G4Material* absorberMat;
  G4Material* gapMat;
  G4Material* defaultMat;

  // Geometry operations  
  G4RotationMatrix *RotXP = new G4RotationMatrix();
  G4RotationMatrix *RotYP = new G4RotationMatrix();
  G4RotationMatrix *RotZP = new G4RotationMatrix();
  G4RotationMatrix *RotXM = new G4RotationMatrix();
  G4RotationMatrix *RotYM = new G4RotationMatrix();
  G4RotationMatrix *RotZM = new G4RotationMatrix();


  G4ThreeVector posBottomAbsorber;
  G4ThreeVector posTopAbsorber;
  G4ThreeVector posLeftRectAbsorber;
  G4ThreeVector posRightRectAbsorber;
  G4ThreeVector posFrontRectAbsorber;
  G4ThreeVector posBackRectAbsorber;
  G4ThreeVector posSourceCan;
  G4ThreeVector posDet1;
  G4ThreeVector posDet2;
 
  // Geometry Parameters
  // brick full length
  G4double brickX;         G4double brickY;         G4double brickZ;
  G4double HalfBrickX;     G4double HalfBrickY;     G4double HalfBrickZ;
  G4double b2bGap;         G4double delta;
  G4double brickYZ;        G4double brickVolume;

  // Lead Array full length
  G4double PbArrayX;       G4double PbArrayY;       G4double PbArrayZ;

  G4double ArrayXBy2;      G4double ArrayYBy2;      G4double ArrayZBy2;
  G4double ArraySegX;      G4double ArraySegY;      G4double ArraySegZ;
  G4double ArraySegX2;     G4double ArraySegY2;     G4double ArraySegZ2;
  G4double PbArrayXY;      G4double PbArrayVolume;  
  G4double worldSizeX;     G4double worldSizeY;     G4double worldSizeZ;

  //----- Here we are using dx, dy and dz as half lengths to minimize computations.
  //---------------------------------------------------------- BoxBottom and BoxTop ------------------------------
  G4double BBDX;    G4double BBDY;    G4double BBDZ;    G4double BTDX;     G4double BTDY;    G4double BTDZ;
  G4double BBCX;    G4double BBCY;    G4double BBCZ;    G4double BTCX;     G4double BTCY;    G4double BTCZ;
  G4double BBBTVol ;
  //--------------------------------------------------------- RectBottom and RectTop ---------------------------
  G4double RLDX;    G4double RLDY;    G4double RLDZ;    G4double RRDX;     G4double RRDY;    G4double RRDZ;
  G4double RLCX;    G4double RLCY;    G4double RLCZ;    G4double RRCX;     G4double RRCY;    G4double RRCZ;
  G4double RLRRVol;
  //---------------------------------------------------------- BoxFront and BoxBack ------------------------------
  G4double BFDX;    G4double BFDY;    G4double BFDZ;    G4double BKDX;     G4double BKDY;    G4double BKDZ;
  G4double BFCX;    G4double BFCY;    G4double BFCZ;    G4double BKCX;     G4double BKCY;    G4double BKCZ;
  G4double BFBKVol;
  //--------------------------------------------------------- BoxCentre -----------------------------------------
  G4double BCDX;    G4double BCDY;    G4double BCDZ;    G4double BCCX;     G4double BCCY;    G4double BCCZ;
  //--------------------------------------------------------- Hole through BoxFront -----------------------------
  G4double PORTD;   G4double PORTL;
  //---------------------------------------------------------- Cylinder -----------------------------------------
  G4double CYLX;    G4double CYLY;   G4double CYLZ;     G4double CYIR;     G4double CYOR;    
  G4double CYLH;    G4double CBH;   


  // Define Volumes for Geometries

  // World
  G4LogicalVolume* worldL;

  // Full Absorber Volume
  G4LogicalVolume* PbArrayL;
  G4VPhysicalVolume* PbArrayP;

  // BottomPlane
  G4VPhysicalVolume* bottomAbsorber;

  // TopPlane
  G4VPhysicalVolume* topAbsorber;

  //LeftPlane
  G4VPhysicalVolume* leftAbsorber;

  //RightPlane
  G4VPhysicalVolume* rightAbsorber;

  //FrontPlane
  G4VPhysicalVolume* frontAbsorber;

  //BackPlane
  G4VPhysicalVolume* backAbsorber;

  //CylindricalSource
  G4VPhysicalVolume* cylSrcCan;

};

  
// inline function
//
inline const G4VPhysicalVolume* DetectorConstruction::GetAbsorberPV() const {
  return fAbsorberPV;
}

inline const G4VPhysicalVolume* DetectorConstruction::GetTargetPV() const {
  return fTargetPV;
}

#endif

