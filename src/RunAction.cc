// $Id: RunAction.cc 2014-11-29 16:07:06Z ABHIJIT $
// //
// /// \file RunAction.cc
// /// \brief Implementation of the RunAction class

#include "B4RunAction.hh"
#include "B4Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


B4RunAction::B4RunAction() : G4UserRunAction()
{ 
    // set printing event number per each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);     

    // Analysis Manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    C4cout << "Using " << analysisManager->GetType() << G4endl;
    // Create directories 
    //analysisManager->SetHistoDirectoryName("histograms");
    //analysisManager->SetNtupleDirectoryName("ntuple");
    analysisManager->SetVerboseLevel(1);
    analysisManager->SetFirstHistoId(1);
     
    // Book histograms, ntuple
     
    // Creating histograms
    analysisManager->CreateH1("1","Edep in absorber", 100, 0., 100*keV);
    analysisManager->CreateH1("2","trackL in absorber", 100, 0., 2*m);
     
    // Creating ntuple
     
    analysisManager->CreateNtuple("Trans-", "Edep and TrackL");
    analysisManager->CreateNtupleDColumn("Eabs");
    analysisManager->CreateNtupleDColumn("Labs");
    analysisManager->FinishNtuple();
}

RunAction::~RUnAction()
{
    delete G4AnalysisManager:Instance();
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
    // Get Analysis Manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

    G4String fileName = "Trans-";
    analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
    G4Analysismanager* analysisManager = G4AnalysisManager::Instance();
    if ( analysisManager->GetH1(1) ) 
    {
        G4cout << "\n ----> print histograms statistic ";
        if(isMaster) 
        {
            G4cout << "for the entire run \n" << G4endl; 
        }
        else 
        {
            G4cout << "for the local thread \n" << G4endl; 
        } 
        
        G4cout << " EAbs : mean = " << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;
    
        G4cout << " EGap : mean = " << G4BestUnit(analysisManager->GetH1(2)->mean(), "Energy") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Energy") << G4endl;
    
        G4cout << " LAbs : mean = " << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;
        G4cout << " LGap : mean = " << G4BestUnit(analysisManager->GetH1(4)->mean(), "Length") 
               << " rms = " << G4BestUnit(analysisManager->GetH1(4)->rms(),  "Length") << G4endl; 
    }

    // save histograms & ntuple
    
    analysisManager->Write();
    analysisManager->CloseFile(); 
}
      
    
     
     

