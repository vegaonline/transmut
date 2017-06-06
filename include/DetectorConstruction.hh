// $Id: DetectorConstruction.hh 2017-05-29 vega $
// // 
// /// \file DetectorConstruction.hh
// /// \brief Definition of the DetectorConstruction class
//
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;
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
  virtual G4VPhysicalVolume* Construct();       // construct geometry of the setup
  //        virtual void ConstructSDandField();
  G4VPhysicalVolume* Construct(); 
  void UpdateGeometry();


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
  auto nPbX = 0;
  auto nPbY = 0;
  auto nPbZ = 0;
  auto nPb = 0;

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
  G4ThreeVector posFrontRectabsorber;
  G4ThreeVector posBackRectAbsorber;
  G4ThreeVector posSourceCan;
  G4ThreeVector posDet1;
  G4ThreeVector posDet2;
 
  // Geometry Parameters
  // brick full length
  G4double brickX;     G4double brickY;     G4double brickZ;
  auto HalfBrickX;     auto HalfBrickY;     auto HalfBrickZ;
  G4double b2bGap;
  auto delta;
  auto brickYZ;        auto brickVolume;

  // Lead Array full length
  G4double PbArrayX;   G4double PbArrayY;   G4double PbArrayZ;

  auto ArrayXBy2;      auto ArrayYBy2;      auto ArrayZBy2;
  auto ArraySegX;      auto ArraySegY;      auto ArraySegZ;
  auto ArraySegX2;     auto ArraySegY2;     auto ArraySegZ2;
  auto PbArrayXY;      auto PbArrayVolume;  
  auto worldSizeX;     auto worldSizeY;     auto worldSizeZ;

  //----- Here we are using dx, dy and dz as half lengths to minimize computations.
  //---------------------------------------------------------- BoxBottom and BoxTop ------------------------------
  auto BBDX;    auto BBDY;    auto BBDZ;    auto BTDX;     auto BTDY;    auto BTDZ;
  auto BBCX;    auto BBCY;    auto BBCZ;    auto BTCX;     auto BTCY;    auto BTCZ;
  auto BBBTVol ;
  //--------------------------------------------------------- RectBottom and RectTop ---------------------------
  auto RLDX;    auto RLDY;    auto RLDZ;    auto RRDX;     auto RRDY;    auto RRDZ;
  auto RLCX;    auto RLCY;    auto RLCZ;    auto RRCX;     auto RRCY;    auto RRCZ;
  auto RLRRVol;
  //---------------------------------------------------------- BoxFront and BoxBack ------------------------------
  auto BFDX;    auto BFDY;    auto BFDZ;    auto BKDX;     auto BKDY;    auto BKDZ;
  auto BFCX;    auto BFCY;    auto BFCZ;    auto BKCX;     auto BKCY;    auto BKCZ;
  auto BFBKVol;
  //--------------------------------------------------------- BoxCentre -----------------------------------------
  auto BCDX;    auto BCDY;    auto BCDZ;    auto BCCX;     auto BCCY;    auto BCCZ;
  //--------------------------------------------------------- Hole through BoxFront -----------------------------
  auto PORTD;   auto PORTL;
  //---------------------------------------------------------- Cylinder -----------------------------------------
  auto CYLX;    auto CYLY;   auto CYLZ;     auto CYIR;     auto CYOR;    auto CYLH;
  auto CBH;   


  // Define Logical Volume for Geometries

  // World
  G4LogicalVolume* worldL;

  // Full Absorber Volume
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

