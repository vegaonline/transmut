// $Id:vegaPB01.cc 2017-05-29 vega $
//
// Abhijit Bhattacharyya, BARC, MUMBAI
//
/// \file vegaPB01.cc
/// \brief Main program 

#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

#include "Randomize.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "FTFP_BERT.hh"


#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " vegaPB01 [-m macro ] [-u UIsession] " << G4endl;
    G4cerr << "   note: -t option is available only for multi-threaded mode."
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  if ( argc > 7 ) {
    PrintUsage();
    return 1;
  }
  
  G4String macro;
  G4String session;

  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else {
      PrintUsage();
      return 1;
    }
  }  
  
  // Choose the Random engine
  //
  //      G4Random::setTheEngine(new CLHEP::RanecuEngine);
  CLHEP::HepRandom::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the default run manager
  //

  G4RunManager * runManager = new G4RunManager;


  // Set mandatory initialization classes
  //
  DetectorConstruction* detConst = new DetectorConstruction();
  runManager->SetUserInitialization(detConst);
  
  //  myPhysicsList* physicsList = new myPhysicsList();
  G4VModularPhysicsList* physicsList = new FTFP_BERT;
  runManager->SetUserInitialization(physicsList);
    
  //  myPrimaryGeneratorAction* myGen = new myPrimaryGeneratorAction();
  //  runManager->SetUserAction(myGen);
  ActionInitialization* actionInit = new ActionInitialization(detConst);
  runManager->SetUserInitialization(actionInit);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {  
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv, session);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    if (ui->IsGUI())
      UImanager->ApplyCommand("/control/execute gui.mac");
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
