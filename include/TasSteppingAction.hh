// 
// 

#ifndef TasSteppingAction_h
#define TasSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

extern G4double Depos[30];
extern G4String detname[200];
extern int detno;
extern int nDetectors;

class TasDetectorConstruction;
class G4Track;

class TasSteppingAction : public G4UserSteppingAction
{
  public:
    TasSteppingAction(TasDetectorConstruction* myDC);
    virtual ~TasSteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

  private:
    TasDetectorConstruction* myDetector;
};

#endif






