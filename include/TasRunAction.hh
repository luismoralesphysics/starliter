//        
// 

#ifndef RunActionTas_h
#define RunActionTas_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
//include libraries for root
#include "TTree.h"
#include "TFile.h"

extern int nDetectors;

extern  TFile *newfile; // = new TFile("geant_out.root","recreate");
extern  TTree * t;// = new TTree("t","output from geant");
extern  TBranch * ebranch;// = t->Branch("etest",&DeposE[0][0],"etest/D");

extern G4double Depos[30], Depos_gauss[30], Depos_tot;
extern int detno, mult;
extern G4double multi, X[10], Y[10], Z[10], distance;
extern G4String detname[200];

class G4Run;

class RunActionTas : public G4UserRunAction
{
  public:
    RunActionTas();
    virtual ~RunActionTas();

  public:
    virtual void BeginOfRunAction(const G4Run* aRun);
    virtual void EndOfRunAction(const G4Run* aRun);

  private:
    G4int runIDcounter;
};

#endif







