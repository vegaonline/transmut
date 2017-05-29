// $Id: RunAction.hh 2014-11-26 14:41:20Z  ABHIJIT $
// // 
// /// \file RunAction.hh
// /// \brief Definition of the RunAction class

#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

/// Run Action Class
// Edepin absorber, Track Length in Absorber

class RunAction : public G4UserRunAction
{
    public:
        RunAction();
        virtual ~RunAction();

        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);
};

#endif
