#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4ios.hh"

//#include "LHEP_PRECO_HP.hh"
#include "QGSP_BIC_HP.hh"

#include "TasDetectorConstruction.hh"
#include "TasPrimaryGeneratorAction.hh"
#include "TasRunAction.hh"
#include "TasEventAction.hh"
#include "TasSteppingAction.hh"

#include <fstream>
#include <iostream>
#include <stdlib.h>

#include "TTree.h"
#include "TFile.h"
#include "TROOT.h"
#include "TRint.h"
#include "TPluginManager.h"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "clover.hh"



int detno;
int mult;
int noOfEvnt;
int noOfEnerg;
int RunTimes;
int jj;
int eve;
float eventArray[8300][8300];
float evtrun[8300];
float gEnergy[10];
G4double Depos[30];
G4double DeposkeV[30];
G4double Depos_gauss[30];
G4double sigma[30];
G4double Depos_tot;
G4double multi;
G4double X[10], Y[10], Z[10];
G4double Xhit[10], Yhit[10], Zhit[10];
G4double distance;
G4String detname[200];

TFile *newfile;
TTree *t;
TBranch *ebranch;

int nDetectors;

int main(int argc, char** argv)
{
  // Construct the default run manager
  G4RunManager *runManager = new G4RunManager;

  // set mandatory initialization classes
  TasDetectorConstruction *TasDet = new TasDetectorConstruction;
  runManager->SetUserInitialization(TasDet);
//  runManager->SetUserInitialization(new LHEP_PRECO_HP);//Library for neutrons
runManager->SetUserInitialization(new QGSP_BIC_HP);


//gROOT->GetPluginManager()->AddHandler("TVirtualStreamerInfo","*","TStreamerInfo","RIO","TStreamerInfo()");

  #ifdef G4VIS_USE
    G4VisManager* visManager = new G4VisExecutive;
     visManager->Initialize();
  #endif

 // set optional user action class
  runManager->SetUserAction(new RunActionTas);


  // Read input file "inpudat" with level schema
  ifstream fin("inputdat");

  G4cout<<"*****************************************************"<<G4endl;

  fin>>RunTimes;
  fin>>noOfEvnt;
  fin>>noOfEnerg;

  int i,j;

  for (i=0; i<noOfEvnt; ++i)
  {
     fin>>eventArray[i][0];
     if (i==0)
       evtrun[i] = eventArray[i][0]*RunTimes/100.0;
     else
       evtrun[i] = evtrun[i-1] + eventArray[i][0]*RunTimes/100.0;


     for (j=1; j<=noOfEnerg; j++)
      {
        fin>>eventArray[i][j];
      }

  }


  // set mandatory user action class
  runManager->SetUserAction(new TasPrimaryGeneratorAction);

  // set optional user action class
  //  runManager->SetUserAction(new RunActionTas);

  TasEventAction *eventAction = new TasEventAction;
  runManager->SetUserAction(eventAction);

  runManager->SetUserAction(new TasSteppingAction(TasDet));


  // Initialize G4 kernel
  runManager->Initialize();

    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (argc!=1) {
        // interactive mode : define UI session
#ifdef G4UI_USE
        G4UIExecutive* ui = new G4UIExecutive(argc, argv,"Qt");
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
#endif
    }


  runManager->BeamOn(RunTimes);

  // job termination

  fin.close();
  newfile->Close();
  G4cout << "### Rootfile closed!" << G4endl;

  #ifdef G4VIS_USE
     delete visManager;
  #endif

  delete runManager;
  G4cout << "### runManager closed!" << G4endl;

  return 0;
}





