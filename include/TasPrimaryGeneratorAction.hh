//
//

#ifndef TasPrimaryGeneratorAction_h
#define TasPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

extern float gEnergy[10];
extern int noOfEnerg, noOfEvnt, RunTimes, eve;
extern float eventArray[8300][8300], evtrun[8300];
extern int jj;

//class TasDetectorConstruction;
class G4ParticleGun;
class G4Event;
//class TasPrimaryGeneratorMessenger;

class TasPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    TasPrimaryGeneratorAction();
    ~TasPrimaryGeneratorAction();

  public:
    void GeneratePrimaries(G4Event *anEvent);
//    void SetRndmFlag(G4String val){rndmFlag = val;};

  private:
    G4ParticleGun* particleGun[10];
  //    G4ParticleGun* particleGun2;
  //  G4ParticleGun* particleGun3;


//    TasPrimaryGeneratorMessenger *gunMessenger;
//    G4String rndmFlag;
};

#endif





