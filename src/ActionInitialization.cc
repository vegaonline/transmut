// $Id: ActionInitialization.cc 2014-11-26 14:47:43Z  ABHIJIT $
//
/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "ActionInitialization.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "DetectorConstruction.hh"


ActionInitialization::ActionInitialization(DetectorConstruction* detConstruction)
                     : G4VUserActionInitialization(),
                     fDetConstruction(detConstruction) { }

ActionInitialization::~ActionInitialization() { }

void ActionInitializtion::BuildForMaster() const
{
    SetUserAction(new RunAction);
}

ActionInitialization::Build() const
{
    SetUserAction(new PrimaryGeneratorAction);
    SetUserAction(new RunAction);
    EventAction* eventAction = new EventAction;
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(fDetConstruction, eventAction));
}


