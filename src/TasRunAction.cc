//
// 

#include "TasRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4String.hh"

#include <stdio.h>
#include <sstream>
#include <iostream>
#include "globals.hh"

G4String parname_ene[200];
G4String parname_eneD[200];
G4String parname_X[200], parname_XD[200];
G4String parname_Y[200], parname_YD[200];
G4String parname_Z[200], parname_ZD[200];


RunActionTas::RunActionTas()
{
  runIDcounter = 0;
}

RunActionTas::~RunActionTas()
{}

void RunActionTas::BeginOfRunAction(const G4Run* aRun)
{
  ((G4Run *)(aRun))->SetRunID(runIDcounter++);
    
    //-------------------
    // Display info about the detector
    
    G4cout << "Number of active crystals: " << nDetectors << G4endl;
    //-------------------

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  srand(time(0));

  //creates the file and tree for the .root output
  newfile = new TFile("/Users/luismorales/Documents/Research/simon/geant4work/rootFiles/output.root","recreate");
  t = new TTree("t","output from geant");

  detno=1;
  while(detno<nDetectors+1)
    {
     parname_ene[detno]=detname[detno]+"_ene";
     parname_eneD[detno]=detname[detno]+"_ene/D";
     //parname_X[detno]=detname[detno]+"_X";
     //parname_XD[detno]=detname[detno]+"_X/D";
     //parname_Y[detno]=detname[detno]+"_Y";
     //parname_YD[detno]=detname[detno]+"_Y/D";
     //parname_Z[detno]=detname[detno]+"_Z";
     //parname_ZD[detno]=detname[detno]+"_Z/D";
  
     ebranch = t->Branch(parname_ene[detno],&Depos_gauss[detno],parname_eneD[detno]);
     //ebranch = t->Branch(parname_X[detno],&X[detno],parname_XD[detno]);
     //ebranch = t->Branch(parname_Y[detno],&Y[detno],parname_YD[detno]);
     //ebranch = t->Branch(parname_Z[detno],&Z[detno],parname_ZD[detno]);
     ++detno;
    }
   ebranch = t->Branch("total_ene",&Depos_tot,"total_ene/D");
   ebranch = t->Branch("multi",&multi,"multi/D");
   //ebranch = t->Branch("dist",&distance,"dist/D");

  if(G4VVisManager::GetConcreteInstance())
   {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis~/clear/view");
     G4UImanager::GetUIpointer()->ApplyCommand("/vis~/draw/current");
   }

}


void RunActionTas::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " ended." << G4endl;

  if(G4VVisManager::GetConcreteInstance())
   {
     G4UImanager::GetUIpointer()->ApplyCommand("/vis~/show/view");
   }

  t->Write(); //write the tree to file


}












