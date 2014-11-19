// 
// 

#ifndef TasEventAction_h
#define TasEventAction_h 1

#include "G4UserEventAction.hh"

#include "globals.hh"

extern int nDetectors;

extern G4double Depos[30], DeposkeV[30], Depos_gauss[30], sigma[30], Depos_tot;
extern int detno, mult, jj, eve, noOfEvnt;
extern G4double multi, X[10], Y[10], Z[10], Xhit[10], Yhit[10], Zhit[10], distance;
extern float evtrun[8300];

class G4Event;

class TasEventAction : public G4UserEventAction
{
  public:
    TasEventAction();
    ~TasEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);

  private:

};

#endif

    






