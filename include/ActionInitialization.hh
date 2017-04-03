// $Id: ActionInitialization.hh 2014-11-26 14:47:43Z ABHIJIT  $
// //
// /// \file ActionInitialization.hh
// /// \brief Definition of the ActionInitialization class
//
#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

class DetectorConstruction;

/// Action initialization class.
///

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


class ActionInitialization : public G4VUserActionInitialization
{
    public:
        ActionInitialization(DetectorConstruction*);
        virtual ~ActionInitialization();

        virtual void BuildForMaster() const;
        virtual void Build() const;
    private:
        DetectorConstruction* fDetConstruction;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

