
#include "TasDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4CutTubs.hh"

#include "G4UnitsTable.hh"
#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <iostream>
#include <sstream>
#include "G4String.hh"
#include "G4ios.hh"
#include <stdio.h>

#include "clover.hh"

TasDetectorConstruction::TasDetectorConstruction()
    : NaI(0), Al(0), N78O21Ar1(0), Cr20Ni8Fe76(0), C2F4(0), Ge(0) {}

TasDetectorConstruction::~TasDetectorConstruction() {}

G4VPhysicalVolume *TasDetectorConstruction::Construct() {
  DefineMaterials();

  return ConstructDetector();
}

void TasDetectorConstruction::DefineMaterials() {

  // define Parameters
  G4String name, symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;

  // define Elements

  a = 22.99 * g / mole;
  G4Element *Na = new G4Element(name = "Sodium", symbol = "Na", z = 11., a);

  a = 126.90 * g / mole;
  G4Element *I = new G4Element(name = "Iodine", symbol = "I", z = 53., a);

  a = 204.38 * g / mole;
  G4Element *Tl = new G4Element(name = "Thalium", symbol = "Tl", z = 81., a);

  a = 26.982 * g / mole;
  G4Element *elAl =
      new G4Element(name = "element_Aluminum", symbol = "elAl", z = 13., a);

  a = 14.00 * g / mole;
  G4Element *N = new G4Element(name = "Nitrogen", symbol = "N", z = 7., a);

  a = 16.00 * g / mole;
  G4Element *O = new G4Element(name = "Oxygen", symbol = "O", z = 8., a);

  a = 39.95 * g / mole;
  G4Element *Ar = new G4Element(name = "Argon", symbol = "Ar", z = 18., a);

  a = 51.996 * g / mole;
  G4Element *Cr = new G4Element(name = "Chromium", symbol = "Cr", z = 24., a);

  a = 58.69 * g / mole;
  G4Element *Ni = new G4Element(name = "Nickel", symbol = "Ni", z = 28., a);

  a = 55.847 * g / mole;
  G4Element *Fe = new G4Element(name = "Iron", symbol = "Fe", z = 26., a);

  a = 12.011 * g / mole;
  G4Element *C = new G4Element(name = "Carbon", symbol = "C", z = 6., a);

  a = 18.998 * g / mole;
  G4Element *F = new G4Element(name = "Fluorine", symbol = "F", z = 9., a);

  a = 72.630 * g / mole;
  G4Element *elGe =
      new G4Element(name = "Germanium", symbol = "elGe", z = 32., a);
  // define Materials

  //..........Stainless Steel..........

  density = 8.0 * g / cm3;
  Cr20Ni8Fe76 =
      new G4Material(name = "Stainless_Steel", density, ncomponents = 3);
  Cr20Ni8Fe76->AddElement(Cr, natoms = 20);
  Cr20Ni8Fe76->AddElement(Fe, natoms = 76);
  Cr20Ni8Fe76->AddElement(Ni, natoms = 8);

  //..........Polytetrafluorine (PTFE)............

  density = 2.20 * g / cm3;
  C2F4 = new G4Material(name = "PTFE", density, ncomponents = 2);
  C2F4->AddElement(C, natoms = 2);
  C2F4->AddElement(F, natoms = 4);

  //..........NaI......................

  density = 3.67 * g / cm3;
  NaI = new G4Material(name = "Sodium Iodine", density, ncomponents = 3);
  NaI->AddElement(Na, natoms = 1000);
  NaI->AddElement(I, natoms = 1000);
  NaI->AddElement(Tl, natoms = 1);

  //..........Al.......................

  density = 2.698 * g / cm3;
  Al = new G4Material(name = "Aluminum", density, ncomponents = 1);
  Al->AddElement(elAl, natoms = 1);

  //..........Ge.......................

  density = 5.323 * g / cm3;
  Ge = new G4Material(name = "Germanium", density, ncomponents = 1);
  Ge->AddElement(elGe, natoms = 1);

  //..........Air.....................

  density = 1.2927 * mg / cm3;
  N78O21Ar1 = new G4Material(name = "Air", density, ncomponents = 3);
  N78O21Ar1->AddElement(N, natoms = 78);
  N78O21Ar1->AddElement(O, natoms = 21);
  N78O21Ar1->AddElement(Ar, natoms = 1);

  // Print out Elements and Materials
  G4cout << "\n\n "
            "####-------------------------------------------------------#### "
            "\n";
  G4cout << "\n\n\n\n\t\t #### List of elements used #### \n";
  G4cout << *(G4Element::GetElementTable());
  G4cout << "\n\n\n\n\t\t #### List of materials used #### \n";
  G4cout << *(G4Material::GetMaterialTable());
  G4cout << "\n\n "
            "####-------------------------------------------------------#### "
            "\n";
}

G4VPhysicalVolume *TasDetectorConstruction::ConstructDetector() {

  //    // Rotation Matries
  G4RotationMatrix *rotZ90 = new G4RotationMatrix();
  rotZ90->rotateZ(90 * deg);
  G4RotationMatrix *rotZ180 = new G4RotationMatrix();
  rotZ180->rotateZ(180 * deg);
  G4RotationMatrix *rotZ270 = new G4RotationMatrix();
  rotZ270->rotateZ(270 * deg);
  G4RotationMatrix *rotX90 = new G4RotationMatrix();
  rotX90->rotateX(90 * deg);

  //------------------ Experimental room
  //--------------------------------------------------------//

  G4Tubs *room_tube = new G4Tubs("room", 0.0 * cm, 100.0 * cm, 300.0 * cm,
                                 0.0 * deg, 360.0 * deg);
  G4LogicalVolume *room_log =
      new G4LogicalVolume(room_tube, N78O21Ar1, "room", 0, 0, 0);
  G4VPhysicalVolume *room_phys =
      new G4PVPlacement(0, G4ThreeVector(0.0 * cm, 0.0 * cm, 0.0 * cm), "room",
                        room_log, NULL, false, 0);

  //------------------ Dimensions of the volumes
  //------------------------------------------------//

  // germanium leaves
  double leafLength = 20 * mm;
  double leafRadius = 15 * mm;
  // ge Casing
  double geCasingLength = 50 * mm;
  double geCasingWidth = 35 * mm;
  double geCasingThickness = 1 * mm;
  double geCasingWindowThickness =
      0.5 * mm; // window thickness (implement later)
  // target chamber
  double chTopThickness = 3.8 * mm;
  double chOuterThickness = 3.8 * mm;
  double chOuterRadius = 580 * mm + chOuterThickness;
  double chHeight = 220 * mm;
  // BGO shield
  double bgoShieldLength = 50 * mm;
  double bgoOuterBoxWidth = 150 * mm; // external width of the shield
  double bgoInnerBoxWidth =
      36 * mm; // opening in the center to hold the Ge detector
  double bgoFrontWidth = 40 * mm; // front of the trapezoid
  double bgoTrdLength = 50 * mm;

  // Target Chamber
  G4double target_chamber_central_cylinder_innner_radius = 50.0 * mm;
  G4double target_chamber_central_cylinder_outer_radius = 52.0 * mm;
  G4double target_chamber_central_cylinder_z_half_length = 50.0 * mm;
  G4double target_chamber_central_cylinder_start_phi_angle = 0 * deg;
  G4double target_chamber_central_cylinder_delta_phi_angle = 360 * deg;

  G4double target_chamber_central_cylinder_beamline_input_innner_radius =
      4.8 * mm;
  G4double target_chamber_central_cylinder_beamline_input_outer_radius =
      5.0 * mm;
  G4double target_chamber_central_cylinder_beamline_input_z_half_length =
      50.0 * mm;
  G4double target_chamber_central_cylinder_beamline_input_start_phi_angle =
      0 * deg;
  G4double target_chamber_central_cylinder_beamline_input_delta_phi_angle =
      360 * deg;

  //--------- Logical volumes
  //------------------------------------------------------------------//

  // leaf
  G4Tubs *geLeaf_tube = new G4Tubs("geLeaf_tube", 0 * mm, leafRadius,
                                   leafLength, 0 * deg, 360 * deg);
  G4LogicalVolume *geLeaf_log =
      new G4LogicalVolume(geLeaf_tube, Ge, "leaf_log", 0, 0, 0);

  // casing
  G4Box *geCasing_box2 =
      new G4Box("geCasing_box2", geCasingWidth, geCasingWidth, geCasingLength);
  G4Box *geCasing_box1 = new G4Box(
      "geCasing_box1", geCasingWidth - geCasingThickness,
      geCasingWidth - geCasingThickness, geCasingLength - geCasingThickness);
  G4SubtractionSolid *geCasing_sub =
      new G4SubtractionSolid("geCasing_sub", geCasing_box2, geCasing_box1);
  G4LogicalVolume *geCasing_log =
      new G4LogicalVolume(geCasing_sub, Al, "geCasing_log", 0, 0, 0);

  // chamber
  G4Tubs *alChamberTop_tube =
      new G4Tubs("alChamberTop_tube", 0 * mm, chOuterRadius, chTopThickness,
                 0 * deg, 360 * deg);
  G4Tubs *alChamberSide_tube =
      new G4Tubs("alChamberSide_tube", chOuterRadius - chOuterThickness,
                 chOuterRadius, chHeight, 0 * deg, 360 * deg);
  G4Tubs *alChamberBottom_tube =
      new G4Tubs("alChamberBottom_tube", 0 * mm, chOuterRadius,
                 chOuterThickness, 0 * deg, 360 * deg);
  // add all the pieces of the chamber together
  G4UnionSolid *alChamber1_union = new G4UnionSolid(
      "alChamber1_union", alChamberSide_tube, alChamberTop_tube, 0,
      G4ThreeVector(0 * mm, 0 * mm, chHeight + 0.5 * chTopThickness));
  G4UnionSolid *alChamber_union = new G4UnionSolid(
      "alChamber_union", alChamber1_union, alChamberBottom_tube, 0,
      G4ThreeVector(0 * mm, 0 * mm, -chHeight - chOuterThickness));
  G4LogicalVolume *alChamber_log =
      new G4LogicalVolume(alChamber_union, Al, "alChamber_log", 0, 0, 0);

  // BGO shield
  G4Box *bgoCube_box2 = new G4Box("bgoCube_box2", bgoOuterBoxWidth,
                                  bgoOuterBoxWidth, bgoShieldLength);
  G4Box *bgoCube_box1 =
      new G4Box("bgoCube_box2", bgoInnerBoxWidth, bgoInnerBoxWidth,
                bgoShieldLength + bgoTrdLength + 1 * mm);
  G4Trd *bgoTrd_trd = new G4Trd("bgoTrd_trd", bgoOuterBoxWidth, bgoFrontWidth,
                                bgoOuterBoxWidth, bgoFrontWidth, bgoTrdLength);
  G4UnionSolid *bgoShield_union = new G4UnionSolid(
      "bgoShield_union", bgoCube_box2, bgoTrd_trd, 0,
      G4ThreeVector(0 * mm, 0 * mm, bgoShieldLength + bgoTrdLength));
  G4SubtractionSolid *bgoShield_sub =
      new G4SubtractionSolid("bgoShield_sub", bgoShield_union, bgoCube_box1, 0,
                             G4ThreeVector(0 * mm, 0 * mm, bgoTrdLength));

  // G4LogicalVolume* bgoCube_log = new
  // G4LogicalVolume(bgoCube_sub,NaI,"bgoCube_log",0,0,0);
  G4LogicalVolume *bgoShield_log =
      new G4LogicalVolume(bgoShield_sub, NaI, "bgoShield_log", 0, 0, 0);

  // Target Chamber
  G4Tubs *target_chamber_central_cylinder =
      new G4Tubs("target_chamber_central_cylinder", 0. * mm,
                 target_chamber_central_cylinder_outer_radius,
                 target_chamber_central_cylinder_z_half_length,
                 target_chamber_central_cylinder_start_phi_angle,
                 target_chamber_central_cylinder_delta_phi_angle);
  // Constructed for the purpose to be used in G4SubtractionSolid.
  G4Tubs *target_chamber_central_cylinder_cut =
      new G4Tubs("target_chamber_central_cylinder_cut", 0. * mm,
                 target_chamber_central_cylinder_innner_radius,
                 target_chamber_central_cylinder_z_half_length + 1 * mm,
                 target_chamber_central_cylinder_start_phi_angle,
                 target_chamber_central_cylinder_delta_phi_angle);
  G4Tubs *target_chamber_central_cylinder_beamline_input = new G4Tubs(
      "target_chamber_central_cylinder_beamline_input",
      target_chamber_central_cylinder_beamline_input_innner_radius,
      target_chamber_central_cylinder_beamline_input_outer_radius,
      target_chamber_central_cylinder_beamline_input_z_half_length,
      target_chamber_central_cylinder_beamline_input_start_phi_angle,
      target_chamber_central_cylinder_beamline_input_delta_phi_angle);
  G4Tubs *target_chamber_central_cylinder_beamline_input_cut = new G4Tubs(
      "target_chamber_central_cylinder_beamline_input_cut", 0,
      target_chamber_central_cylinder_beamline_input_innner_radius,
      target_chamber_central_cylinder_beamline_input_z_half_length + 1,
      target_chamber_central_cylinder_beamline_input_start_phi_angle,
      target_chamber_central_cylinder_beamline_input_delta_phi_angle);
  G4LogicalVolume *target_chamber_central_cylinder_logical =
      new G4LogicalVolume(target_chamber_central_cylinder, Al,
                          "target_chamber_central_cylinder", 0, 0, 0);
  G4LogicalVolume *target_chamber_central_cylinder_beamline_input_logical =
      new G4LogicalVolume(target_chamber_central_cylinder_beamline_input, Al,
                          "target_chamber_central_cylinder_beamline_input", 0,
                          0, 0);
  G4VSolid *target_chamber_cylinder_union1 = new G4UnionSolid(
      "target_chamber_cylinder_union1", target_chamber_central_cylinder,
      target_chamber_central_cylinder_beamline_input, rotX90,
      G4ThreeVector(0. * mm, 0. * mm, 0. * mm));
  G4VSolid *target_chamber_cylinder_union2 = new G4SubtractionSolid(
      "target_chamber_cylinder_union2", target_chamber_cylinder_union1,
      target_chamber_central_cylinder_cut);
  G4VSolid *target_chamber_cylinder_union = new G4SubtractionSolid(
      "target_chamber_cylinder_union", target_chamber_cylinder_union2,
      target_chamber_central_cylinder_beamline_input_cut, rotX90,
      G4ThreeVector(0. * mm, 0. * mm, 0. * mm));

  G4LogicalVolume *target_chamber_cylinder_union_logical =
      new G4LogicalVolume(target_chamber_cylinder_union, Al,
                          "target_chamber_cylinder_union", 0, 0, 0);

  //--------- Physical volumes
  //------------------------------------------------------------------//

  G4RotationMatrix *rotY90 = new G4RotationMatrix();
  rotY90->rotateY(90 * deg);

  G4RotationMatrix *rotY90_rotZ90 =
      new G4RotationMatrix(0 * deg, 90 * deg, 0 * deg);
  // coordinates of the center of the clover window
  double posX = 470 * mm;
  double posY = 0 * mm;
  double posZ = 0 * mm;

  // bgo shield offset
  double shieldX = chOuterRadius;

  // leaf offset within clover (distance in X and Y from the clover center)
  double leafX = leafRadius;
  double leafY = leafRadius;

  detno = 0;
  std::ostringstream oss;
  G4String s;

  // chamber
  new G4PVPlacement(rotX90, G4ThreeVector(0 * mm, 0 * mm, 0 * mm),
                    alChamber_log, "alChamber", room_log, false, 0);

  // clover casing
  new G4PVPlacement(rotY90, G4ThreeVector(posX + geCasingLength, posY, posZ),
                    geCasing_log, "geCasing", room_log, false, 0);

  // BGO shield
  new G4PVPlacement(rotY90, G4ThreeVector(posX + shieldX, posY, posZ),
                    bgoShield_log, "bgoShield", room_log, false, 0);
  // G4VPhysicalVolume* bgoTrd_phys = new
  // G4PVPlacement(0,G4ThreeVector(0*mm,0*mm,bgoShieldLength+bgoTrdLength),bgoTrd_log,"bgoTrd",room_log,false,0);

  // Target Chamber
  new G4PVPlacement(rotY90_rotZ90, G4ThreeVector(0 * mm, 0 * mm, 0 * mm),
                    target_chamber_cylinder_union_logical,
                    "target_chamber_cylinder_union", room_log, false, 0);

  detno = 0;
  while (detno < 5) {
    detno++;
    oss.str(""); // clear the ostringstream
    oss << detno; // stream int into ostringstream
    s = oss.str(); // copy ostringstream to G4string
    detname[detno] = "NaI" + s;
  }

  new G4PVPlacement(rotY90, G4ThreeVector(posX + geCasingLength, leafY, leafX),
                    geLeaf_log, detname[1], room_log, false, 0);
  new G4PVPlacement(rotY90, G4ThreeVector(posX + geCasingLength, leafY, -leafX),
                    geLeaf_log, detname[2], room_log, false, 0);
  new G4PVPlacement(rotY90, G4ThreeVector(posX + geCasingLength, -leafY, leafX),
                    geLeaf_log, detname[3], room_log, false, 0);
  new G4PVPlacement(rotY90,
                    G4ThreeVector(posX + geCasingLength, -leafY, -leafX),
                    geLeaf_log, detname[4], room_log, false, 0);
  /*
          detno++;
          oss.str("");   //clear the ostringstream
          oss << detno;  //stream int into ostringstream
          s=oss.str();   //copy ostringstream to G4string
          detname[detno]="NaI" + s;

          // ----- NaI cube
          G4VPhysicalVolume* crystal_phys1 = new
     G4PVPlacement(0,G4ThreeVector(posX,posY,posZ),crystal_log,detname[detno],room_log,false,0);


  */

  // total number of detectors (to make sure all is saved in ROOT tree and
  // summed properly
  detno = 4;
  nDetectors = detno;

//   Clover *aDetector = new Clover();
//   aDetector->SetPosition(G4ThreeVector(0. * mm, 0. * mm, 0. * mm));
//   G4RotationMatrix rotMat;
//   rotMat.set(0, 0, 0);
//   rotMat.invert();
//   aDetector->SetRotation(rotMat);
//   aDetector->Placement(0, room_phys);

  //========================== Visualization attributes
  //=========================================//

  room_log->SetVisAttributes(G4VisAttributes::Invisible);

  // visualization for scintillators = GREEN
  G4VisAttributes *GreenAttr = new G4VisAttributes(G4Colour(0., 1., 0.));
  GreenAttr->SetVisibility(true);
  GreenAttr->SetForceSolid(true);

  // visualization for reflector = PURPLE
  G4VisAttributes *PurpleAttr = new G4VisAttributes(G4Colour(1., 0., 1.));
  PurpleAttr->SetVisibility(true);
  PurpleAttr->SetForceSolid(true);

  // visualization for aluminum = GREY
  G4VisAttributes *GreyAttr = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
  GreyAttr->SetVisibility(true);
  GreyAttr->SetForceSolid(false);

  // visualization for BLUE
  G4VisAttributes *BlueAttr = new G4VisAttributes(G4Colour(0., 0., 1.));
  BlueAttr->SetVisibility(true);
  BlueAttr->SetForceSolid(true);

  // visualization for RED
  G4VisAttributes *RedAttr = new G4VisAttributes(G4Colour(1., 0., 0.));
  RedAttr->SetVisibility(true);
  RedAttr->SetForceSolid(true);

  // applying the color scheme to the logical volumes

  bgoShield_log->SetVisAttributes(GreenAttr);
  geLeaf_log->SetVisAttributes(PurpleAttr);
  geCasing_log->SetVisAttributes(BlueAttr);
    alChamber_log->SetVisAttributes(GreyAttr);
    target_chamber_cylinder_union_logical->SetVisAttributes(GreenAttr);

 /*
    inset_log->SetVisAttributes(GreenAttr);
    inset_reflective_log->SetVisAttributes(PurpleAttr);
    inset_al_log->SetVisAttributes(GreyAttr);
 */
//always return the world at the end
  return room_phys;
}


