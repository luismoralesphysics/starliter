//
#ifndef TasDetectorConstruction_h
#define TasDetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
  
class G4VPhysicalVolume;
class G4Material;

extern  G4String detname[200];
extern int detno;
extern G4double X[10], Y[10], Z[10];

extern int nDetectors;


class TasDetectorConstruction : public G4VUserDetectorConstruction{
   public:
      TasDetectorConstruction();
      ~TasDetectorConstruction();

   public:
    G4VPhysicalVolume* Construct();

   private:
      void DefineMaterials();
  G4Material *NaI, *Al, *N78O21Ar1, *Cr20Ni8Fe76, *C2F4, *Ge;

      G4VPhysicalVolume* ConstructDetector();

};
 

#endif

