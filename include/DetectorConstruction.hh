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
        virtual G4VPhysicalVolume* Construct();
  //        virtual void ConstructSDandField();

        //get Methods
        //
      const G4VPhysicalVolume* GetAbsorberPV() const;       // For Lead Blocks to epi-thermalize
         const G4VPhysicalVolume* GetTargetPV() const;         // For Be target to produce neutron

    private:
        // Methods
        //
        void DefineMaterials();
        G4VPhysicalVolume* DefineVolumes();

        // Data members
        //
        // static G4ThreadLocal G4GlobalMagFieldMessenger* fMagFieldMessenger;
        G4VPhysicalVolume* fAbsorberPV;           // the absorber physical volume
        G4VPhysicalVolume* fTargetPV;             // the Be target to produce neutron

        G4bool fCheckOverlaps;                   // check volumes overlap
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

