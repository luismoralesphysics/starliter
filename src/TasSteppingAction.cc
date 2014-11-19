// 
// 
#include "TasSteppingAction.hh"

#include "TasDetectorConstruction.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4ios.hh"
#include <stdio.h>
#include <cstring>

TasSteppingAction::TasSteppingAction(TasDetectorConstruction* myDC):myDetector(myDC)
{ }

void TasSteppingAction::UserSteppingAction(const G4Step* aStep)
{

 const G4VPhysicalVolume* currentVolume1 = aStep->GetPreStepPoint()-> GetPhysicalVolume();

 if (currentVolume1 != NULL)
 {  
   detno=1;
   while (detno<nDetectors+1)
    {
      if (currentVolume1->GetName() == detname[detno])   
      {
        Depos[detno] += aStep->GetTotalEnergyDeposit();
      }

      ++detno;
    }
 }

}









