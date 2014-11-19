//    
//

#include "TasPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh" 
#include "G4ParticleDefinition.hh"      

#include "G4UImanager.hh"

#include "globals.hh"
#include "Randomize.hh"

#include <iostream>
#include <fstream>
#include <iomanip>


TasPrimaryGeneratorAction::TasPrimaryGeneratorAction()
{
  G4int n_particle = 1;

  // Setting the default particle
  G4ParticleTable *pTable = G4ParticleTable::GetParticleTable();
  G4String pName;
  G4ParticleDefinition *particle = pTable->FindParticle(pName="gamma");

  int gno1=1;

  for(gno1=0; gno1<=noOfEnerg; ++gno1)
    {
       particleGun[gno1] = new G4ParticleGun(n_particle);

       particleGun[gno1]->SetParticleDefinition(particle);
    }

}


TasPrimaryGeneratorAction::~TasPrimaryGeneratorAction()
{
  int gno1;
  for(gno1=0 ; gno1<noOfEnerg; ++ gno1)
    {
  delete particleGun[gno1];
 }
}

void TasPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)    
{
 G4double px[15], py[15], pz[15], costheta[15], sintheta[15], phi[15];       // CoM frame
 G4double pxLab[15], pyLab[15], pzLab[15], cosLabTheta[15], sinLabTheta[15]; // Lab frame
 G4double dopplerEffect[15], gmma[15], velocity[15];                         // doppler variables
 G4double beamE, mamu;                                                       // energy variables
 int gno=0; 
 G4double x[15], y[15], z[15];

  mamu= 931.5 *MeV;
  //beamE=0.9834; //*MeV/u
  beamE=0.0; //*MeV/u

  for(gno=1;gno<=noOfEnerg; ++gno)
   {

    px[gno] = 0.0;
    py[gno] = 0.0;
    pz[gno] = 0.0;

    costheta[gno]=2.0*G4UniformRand()-1.0;                           // generates a random value for costheta between -1 and 1
    sintheta[gno]=std::sqrt(1.0-costheta[gno]*costheta[gno]);        // sin^2 + cos^2 = 1
    phi[gno]=twopi* G4UniformRand();                                 // generates a random angle between 0 and 2pi
    px[gno]=sintheta[gno] * std::cos(phi[gno]);                      // calculates x = sintheta * cosphi
    py[gno]=sintheta[gno] * std::sin(phi[gno]);                      // calculates y = sintheta * sinphi
    pz[gno]=costheta[gno];                                           // calculates z = costheta

    // Now we start necessary calcs for Doppler effects
     gmma[gno] = ((beamE*MeV)/mamu) + 1.;                            // E_kinetic = gamma*mass*c*c - mass*c*c
     velocity[gno] = sqrt(1. - (1./(gmma[gno]*gmma[gno])));          // gamma = 1/sqrt(1-(v/c)^2)


   // Move to the lab frame (we assume here that the beam axis is the z axis)
     cosLabTheta[gno] = (costheta[gno] + velocity[gno])/(1. + velocity[gno]*costheta[gno]);   // relativistic effects
     sinLabTheta[gno]=std::sqrt(1.0-cosLabTheta[gno]*cosLabTheta[gno]);                       // sin^2 + cos^2 = 1
     pxLab[gno]=sinLabTheta[gno] * std::cos(phi[gno]);                                        // calculates x = sintheta * cosphi
     pyLab[gno]=sinLabTheta[gno] * std::sin(phi[gno]);                                        // calculates y = sintheta * sinphi
     pzLab[gno]=cosLabTheta[gno];                                                             // calculates z = costheta
     dopplerEffect[gno] = 1./(gmma[gno]*(1. - velocity[gno]*cosLabTheta[gno]));                // f = 1/gamma * 1/(1-v/c*costheta) * f0


   // Sets the particles momentum with the corrections
     particleGun[gno]->SetParticleMomentumDirection(G4ThreeVector(pxLab[gno], pyLab[gno], pzLab[gno]));

   // Sets the original postion of the particle
     x[gno]=0.0*cm;
     y[gno]=0.0*cm;
     z[gno]=0.0*cm;

     particleGun[gno]->SetParticlePosition(G4ThreeVector(x[gno], y[gno], z[gno]));
 
  // Sets particle energies with the doppler boost
     particleGun[gno]->SetParticleEnergy(dopplerEffect[gno]*eventArray[eve][gno]*keV);

  // To generate an event
     particleGun[gno]->GeneratePrimaryVertex(anEvent);

   }

}















