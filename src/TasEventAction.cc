//
//

#include "TasEventAction.hh"
#include "TasRunAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include <stdio.h>
#include "Randomize.hh"

int hit, i=0;

TasEventAction::TasEventAction() {}

TasEventAction::~TasEventAction() {}

void TasEventAction::BeginOfEventAction(const G4Event *evt) {
  Depos_tot = 0.0 * MeV;
  mult = 0;
  detno = 1;

  while (detno < nDetectors + 1) // set energies back to zero
  {
    Depos[detno] = 0.0 * eV;
    sigma[detno] = 0.0 * keV;
    DeposkeV[detno] = 0.0 * keV;
    Depos_gauss[detno] = 0.0 * keV;
    ++detno;
  }
}

void TasEventAction::EndOfEventAction(const G4Event* evt)
{
  if ((evt->GetEventID()+1) % 10000 == 0) // prints out event number to screen
   G4cout << ">>> Event " << evt->GetEventID()+1 << G4endl;


 // figures out what gamma cascade to run
 // if the EventID number is large enough, go onto next cascade
  jj=evt->GetEventID()+1;
  for (i=0; i<noOfEvnt; ++i)
  {
     if (i==0 && jj<=evtrun[0])
       eve=i;
     else if (i>0 && jj<=evtrun[i] && jj>evtrun[i-1])
       eve=i;
  }


  detno=1;
  hit=1;


  while (detno<nDetectors+1)
  {
      if (Depos[detno]>0.0)
      {
          DeposkeV[detno]=1000*Depos[detno];

          sigma[detno] = 1;
          //sigma[detno]=-5.59375e-15*pow(DeposkeV[detno],4.0)+1.85975e-10*pow(DeposkeV[detno],3.0)-2.47836e-6*pow(DeposkeV[detno],2.0)+2.33408e-2*DeposkeV[detno]+7.00328;

          Depos_gauss[detno] = G4RandGauss::shoot(DeposkeV[detno],sigma[detno]);

          //if (Depos_gauss[detno] < 160.0) Depos_gauss[detno]=0.0;
          if (Depos_gauss[detno] < 0.0)
          {
              Depos_gauss[detno]=0.0;
          }
          else
          {
              Depos_tot += Depos_gauss[detno];
              ++mult;
              ++hit;
          }
//       G4cout << detno << " " << DeposkeV[detno] << " " << Depos_gauss[detno] << G4endl;
      }
      ++detno;
  }


    multi=mult; // multi is the number of detectors hit

    t->Fill();
}


