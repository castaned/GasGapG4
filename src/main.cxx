#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "G4UImanager.hh"
#include "PhysicsList.hh"
#include "FTFP_BERT_HP.hh"

//#include "G4VModularPhysicsList.hh"

#include "DetectorConstruction.h"
//#include "ActionInitialization.h"
#include "PrimaryGenerator.h"
#include "GasGapRunAction.hh"
#include "GasGapEventAction.hh"
#include "GasGapSteppingAction.hh"
#include "GasGapElectricFieldSetup.hh"
//#include "HistoManager.hh"

#include "G4ScoringManager.hh"

int main(int argc, char** argv)
{

// #ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
// #else
//   G4RunManager* runManager = new G4RunManager;
// #endif

  G4RunManager * runManager = new G4RunManager;
  
  // Construct the helper class to manage the electric field &                                                         // the parameters for the propagation of particles in it.                                                           
  //  GasGapElectricFieldSetup* field = new GasGapElectricFieldSetup() ;    
  
  //Initialize detector construction
  G4VUserPhysicsList* physics = new FTFP_BERT_HP();
  runManager->SetUserInitialization(physics);
  //  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserAction(new PrimaryGenerator) ;
  runManager->SetUserInitialization(new DetectorConstruction());
  
	
  //Initialize physics list
  //	runManager->SetUserInitialization(new PhysicsList());

  
  //Initialize kernel
  //  runManager->Initialize();
	
  //Initialize user action
  GasGapRunAction* run_action = new GasGapRunAction();
  runManager->SetUserAction(run_action) ;         
  runManager->SetUserAction(new GasGapEventAction()) ;
  runManager->SetUserAction(new GasGapSteppingAction);
  //	runManager->SetUserInitialization(new ActionInitialization());
  //	GasGapRunAction* run_action = new GasGapRunAction(histo);
  //	F02RunAction* run_action = new F02RunAction();
  //	F02RunAction* run_action = new F02RunAction(histo);
  //	runManager->SetUserAction(new GasGapEventAction(run_action,histo)) ;
	
  //Activate command-based scorer
  G4ScoringManager::GetScoringManager();
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute share/init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute share/init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  //  delete field;
  
  return 0;
}
