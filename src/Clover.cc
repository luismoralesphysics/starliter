//---------------------------------------------------------------------
// Create the solids defining an Eurogam Phase-II detector
//---------------------------------------------------------------------
#include "Clover.hh"
//#include "MyMaterials.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "globals.hh"
#include "G4NistManager.hh"

#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

Clover::Clover() {
  // some default Clover detector parameters
  fTotalGeL = 69.00 * mm; // was 70
  fCrystalR = 24.25 * mm; // was 25
  fEndCap2Ge = 20.00 * mm; // Distance from Al outer face to Ge
  // added to fudge PTG's efficiency measurements for STUK
  fFudge = 5.0 * mm;
  fEndCap2Ge += fFudge;

  fHoleR = 5.5 * mm; // was 5.0
  fPassiveThick = 0.5 * mm;
  fContactThick = 0.5 * mm;
  fGapBetweenLeaves = 0.6 * mm;

  //---------------------------------
  // make the required materials and assign them
  // fMat = MyMaterials::GetInstance();
  G4NistManager *man = G4NistManager::Instance();
  man->SetVerbose(1);

  // assign default materials.....
  geMaterial =
      man->FindOrBuildMaterial("G4_Ge"); // fMat->GetMaterial("Germanium");
  endCapMaterial =
      man->FindOrBuildMaterial("G4_Al"); // fMat->GetMaterial("DurAl");
  vacuumMaterial =
      man->FindOrBuildMaterial("G4_Galactic"); // fMat->GetMaterial("Vacuum");
  contactMaterial =
      man->FindOrBuildMaterial("G4_Li"); // fMat->GetMaterial("Lithium");

  // create the solids.....
  CreateSolids();
}

// Destructor
Clover::~Clover() {
  //
}

void Clover::SetPosition(G4ThreeVector thisPos) {

  position = thisPos * mm;
  G4cout << " ----> A Phase-II will be placed at " << position / mm << " mm"
         << G4endl;
}

void Clover::SetRotation(G4RotationMatrix thisRot) {

  rotation = thisRot;

  // G4cout << " ----> Single Ge will be placed at " << position/mm << " mm" <<
  // G4endl;
}

//---------------------------------------------------------------------
// Create the solids defining Phase-II Clovers
//---------------------------------------------------------------------
void Clover::CreateSolids() {
  // An approximate CloverII
  G4cout << "Constructing archetypal Clover Clover" << G4endl;

  //---------------------------------------------------------
  // end-cap
  G4double endCapFrontThickness =
      1.5 * mm; // 1.5 in the Duchene paper //1.2*mm; in another simulation
  G4double endCapTaperThickness = 1.5 * mm;
  G4double endCapSideThickness = 1.5 * mm;

  G4double GeGap = fEndCap2Ge; // outsideface of endcap to Ge
  G4double taperAngle = 7.0 * degree;

  G4double endCapTotalL =
      fTotalGeL + GeGap + endCapFrontThickness + 5. * mm; //+ Gap at rear end
  G4double endCapFrontD = 43.5 * mm;
  G4double endCapBackD = 50.5 * mm;
  G4double endCapTaperL = 55.0 * mm;

  G4double endCapBoxL = endCapTotalL - endCapTaperL;

  // the tapered part
  G4Trap *solidTaperedCloverEC = new G4Trap(
      "taperedCloverEC",
      endCapTaperL / 2., // Half z-length [pDz]
      0.00 *
          degree, // Polar angle of line joining centres of the faces @ -/+pDz
      2. * taperAngle, // aequivalent zimuthal angle
      endCapFrontD, // pDy1 half y length at -pDz
      endCapFrontD, // pDx1 half x length at -pDz, -pDy1
      endCapFrontD, // pDx2 half x length at -pDz, +pDy1
      0.00 * degree, // pAlpha1 wrt y-axis from the centre of the side (lower
                     // endcap)
      endCapBackD, // pDy2 half y length at +pDz
      endCapBackD, // pDx3 half x length at +pDz, -pDy2
      endCapBackD, // pDx4 half x length at +pDz, +pDy2
      0.00 * degree); // pAlpha2 wrt y-axis from the centre of the side (upper
                      // endcap)

  // the rectangular part.....
  G4Box *endCapBox =
      new G4Box("endCapBox", endCapBackD, endCapBackD, endCapBoxL / 2.);
  G4ThreeVector transECBox(0. * mm, 0. * mm,
                           endCapTaperL / 2. + endCapBoxL / 2.);

  // add the two together
  solidEndCap = new G4UnionSolid("Box+Taper", solidTaperedCloverEC, endCapBox,
                                 0, transECBox);
  // need the taperL for placement
  fEndCapTaperL = endCapTaperL;

  //---------------------------------------------------------
  // end-cap inner vacuum
  G4double endCapDelta_1 = endCapTaperThickness / cos(taperAngle) -
                           endCapFrontThickness * tan(taperAngle);
  G4double endCapDelta_2 =
      (endCapSideThickness -
       (endCapTaperThickness * sin(taperAngle) * tan(taperAngle) +
        endCapTaperThickness * cos(taperAngle))) /
      tan(taperAngle);

  G4cout << endCapDelta_1 << " " << endCapDelta_2 << G4endl;

  G4double endCapVacTaperL =
      endCapTaperL - endCapFrontThickness; // - endCapDelta_2;
  G4double endCapVacBoxL = endCapBoxL - endCapFrontThickness;
  G4double endCapVacTotalL = endCapVacBoxL + endCapVacTaperL;
  G4double endCapVacFrontD = endCapFrontD - endCapDelta_1;
  G4double endCapVacBackD = endCapBackD - endCapSideThickness;

  // position of vacuum wrt end-cap
  fVacuumPosZ =
      (-endCapTotalL + endCapVacTotalL) / 2. + 1.5 * endCapFrontThickness;

  // tapered part...
  G4Trap *solidTaperVac = new G4Trap(
      "cloverTaperVac",
      endCapVacTaperL / 2., // Half z-length [pDz]
      0.00 *
          degree, // Polar angle of line joining centres of the faces @ -/+pDz
      2. * taperAngle, // aequivalent zimuthal angle
      endCapVacFrontD, // pDy1 half y length at -pDz
      endCapVacFrontD, // pDx1 half x length at -pDz, -pDy1
      endCapVacFrontD, // pDx2 half x length at -pDz, +pDy1
      0.00 * degree, // pAlpha1 wrt y-axis from the centre of the side (lower
                     // endcap)
      endCapVacBackD, // pDy2 half y length at +pDz
      endCapVacBackD, // pDx3 half x length at +pDz, -pDy2
      endCapVacBackD, // pDx4 half x length at +pDz, +pDy2
      0.00 * degree); // pAlpha2 wrt y-axis from the centre of the side (upper
                      // endcap)

  // rectangular part
  G4Box *endCapVacBox = new G4Box("endCapBox", endCapVacBackD, endCapVacBackD,
                                  endCapVacBoxL / 2.);
  G4ThreeVector transVacBox(
      0. * mm, 0. * mm,
      (endCapVacTaperL / 2. + endCapVacBoxL / 2. - 0.0001 * mm));

  // add them together
  solidVacuum = new G4UnionSolid("Vac_Box+Taper", solidTaperVac, endCapVacBox,
                                 0, transVacBox);

  //---------------------------------------------------------
  // The Ge crystal...
  G4double GeTaperL = 36.0 * mm;
  G4double GeTotalL = fTotalGeL; // 70.0 * mm;
  G4double smallSquare = 41.0 * mm;
  G4double largeSquare = 45.5 * mm;

  G4double transX = (largeSquare - smallSquare) / 2.;
  G4double transY = (largeSquare - smallSquare) / 2.;
  fHole_dX = transX + 5. * mm; // 5 mm is a fudge
  fHole_dY = transY + 5. * mm; // 5 mm is a fudge

  // tapered part......
  G4Trap *solidTaper = new G4Trap(
      "cloverTaper",
      GeTaperL / 2., // Half ? z-length [pDz]
      5.05 *
          degree, // Polar angle of line joining centres of the faces @ -/+pDz
      45. * degree, // equivalent azimuthal angle  //DOES NOT MAKE SENSE !!
      smallSquare / 2., // pDy1 half y length at -pDz
      smallSquare / 2., // pDx1 half x length at -pDz, -pDy1
      smallSquare / 2., // pDx2 half x length at -pDz, +pDy1
      0.00 * degree, // pAlpha1 wrt y-axis from the centre of the side (lower
                     // endcap)
      largeSquare / 2., // pDy2 half y length at +pDz
      largeSquare / 2., // pDx3 half x length at +pDz, -pDy2
      largeSquare / 2., // pDx4 half x length at +pDz, +pDy2
      0.0 * degree); // pAlpha2 wrt y-axis from the centre of the side (upper
                     // endcap)

  // HERE !!
  const G4int numZPlanesGe = 4; // no. polycone planes

  G4double zPlaneGe[numZPlanesGe]; // positions of planes
  zPlaneGe[0] = 0.00 * mm;
  zPlaneGe[1] = 2.46 * mm;
  zPlaneGe[2] = 5.00 * mm;
  zPlaneGe[3] = GeTaperL;

  G4double rInnerGe[numZPlanesGe]; // interior radii
  rInnerGe[0] = rInnerGe[1] = rInnerGe[2] = rInnerGe[3] = 0.0 * mm;
  G4double rOuterGe[numZPlanesGe]; // exterior radii
  rOuterGe[0] = 20.5 * mm;
  rOuterGe[1] = 23.54 * mm;
  rOuterGe[2] = rOuterGe[3] = fCrystalR;

  G4Polycone *solidCone =
      new G4Polycone("cloverCone", 0.0 * degree, 360.0 * degree, numZPlanesGe,
                     zPlaneGe, rInnerGe, rOuterGe);

  G4ThreeVector transGeCone(-transX / 2., -transY / 2., -GeTaperL / 2.);
  G4IntersectionSolid *taperedCone = new G4IntersectionSolid(
      "Taper+Cone", solidTaper, solidCone, 0, transGeCone);

  // back part....
  G4double geBoxL = fTotalGeL - GeTaperL;

  G4Box *GeBox =
      new G4Box("GeBox", largeSquare / 2., largeSquare / 2., geBoxL / 2.);
  G4Tubs *GeCyl = new G4Tubs("GeCyl", 0.0 * mm, fCrystalR, geBoxL / 2.,
                             0. * degree, 360. * degree);

  G4ThreeVector transGeBox(transX, transY, 0.0 * mm);
  G4IntersectionSolid *backPart =
      new G4IntersectionSolid("Box+Cyl", GeCyl, GeBox, 0, transGeBox);

  // add front and back
  G4ThreeVector transBack(-transX / 2., -transY / 2.,
                          (GeTaperL / 2. + geBoxL / 2.));
  solidGeLeaf =
      new G4UnionSolid("germanium", taperedCone, backPart, 0, transBack);

  // z-position of Ge-leaf wrt vacuum
  fGeLeafPosZ =
      -endCapVacTaperL / 2. + GeTaperL / 2. + GeGap - endCapFrontThickness;

  G4cout << "end-cap : box/2 " << endCapBoxL / 2. << " taper/2 "
         << endCapTaperL / 2. << " total/2 " << endCapTotalL << G4endl;
  G4cout << "vacuum  : box/2 " << endCapVacBoxL / 2. << " taper/2 "
         << endCapVacTaperL / 2. << " total/2 " << endCapVacTotalL << G4endl;
  G4cout << "ge      : box/2 " << geBoxL / 2. << " taper/2 " << GeTaperL / 2.
         << " total/2 " << GeTotalL << G4endl;

  //------------------------------------------------------------------
  // Inner bore hole + lithium contact + passivated Ge
  G4double GeDepth = 15.00 * mm; // Hole dirilled to this far from face
  G4double passiveThick = fPassiveThick; // passivated Ge
  G4double contactThick = fContactThick; // Li contact

  G4double innerRHole = 0.00 * mm;
  G4double holeR = fHoleR;
  G4double contactR = holeR + contactThick;
  G4double passiveR = contactR + passiveThick;
  G4double holeL = fTotalGeL - GeDepth;
  G4double tubeL = holeL - holeR;

  // the same translation works for all the following rounded tubes
  G4ThreeVector transSphere(
      0.01 * mm, 0.01 * mm,
      -tubeL / 2. - 0.1 * mm); // if offsets are 0. it does not display !!

  // now add a passivated layer
  G4Sphere *passivatedSphere =
      new G4Sphere("passSphere", 0.0 * mm, passiveR, 0. * degree, 360. * degree,
                   0. * degree, 180. * degree);
  G4Tubs *passivatedTube = new G4Tubs("passTube", 0.0 * mm, passiveR,
                                      tubeL / 2., 0. * degree, 360. * degree);
  solidPassivated = new G4UnionSolid("passivatedGe", passivatedTube,
                                     passivatedSphere, 0, transSphere);

  // and the Li contact
  G4Sphere *contactSphere =
      new G4Sphere("sphere1", 0.0 * mm, contactR, 0. * deg, 360. * deg,
                   0. * deg, 180. * deg);
  G4Tubs *contactTube =
      new G4Tubs("tube1", 0.0 * mm, contactR, tubeL / 2., 0. * deg, 360. * deg);
  solidContact =
      new G4UnionSolid("liContact", contactTube, contactSphere, 0, transSphere);

  // bore out a hole
  G4Sphere *boreSphere =
      new G4Sphere("boreSphere", 0.0 * mm, holeR, 0. * degree, 360. * degree,
                   0. * degree, 180. * degree);
  G4Tubs *boreTube = new G4Tubs("boreTube", 0.0 * mm, holeR, tubeL / 2.,
                                0. * degree, 360. * degree);
  solidBoreHole =
      new G4UnionSolid("boreHole", boreTube, boreSphere, 0, transSphere);

  // save this for placement
  fContact_dZ = holeL / 2. - contactThick; // - passiveThick;

  // put corners @ (0,0)
  fGeLeaf_dX = largeSquare / 2. - transX / 2.;
}

//------------------------------------------------------------------
void Clover::Placement(G4int copyNo, G4VPhysicalVolume *physiMother) {
  G4double vacuum_PosZ = fVacuumPosZ;
  G4double geLeaf_PosZ = fGeLeafPosZ;

  logicEndCap =
      new G4LogicalVolume(solidEndCap, endCapMaterial, "clover_EC", 0, 0, 0);
  logicVacuum =
      new G4LogicalVolume(solidVacuum, vacuumMaterial, "clover_Vac", 0, 0, 0);

  for (G4int l = 0; l < 4; l++) {
    logicGeLeaf[l] =
        new G4LogicalVolume(solidGeLeaf, geMaterial, "clover_Leaf", 0, 0, 0);
    logicPassivated[l] =
        new G4LogicalVolume(solidContact, geMaterial, "passivatedGe", 0, 0, 0);
    logicContact[l] = new G4LogicalVolume(solidPassivated, contactMaterial,
                                          "inner_contact", 0, 0, 0);
    logicBoreHole[l] = new G4LogicalVolume(solidBoreHole, vacuumMaterial,
                                           "bore_hole", 0, 0, 0);
  }
  //=================================================================================
  // Do not know why, but the positioning seems to be with respect to the
  // Taper-part :
  // setting the z-position as endCapTaperL/2 puts the front face at z = 0 mm
  //=================================================================================
  // position.setZ(position.z() + fEndCapTaperL/2);

  physiEndCap = new G4PVPlacement(
      &rotation, G4ThreeVector(position.x(), position.y(),
                               position.z() + fEndCapTaperL / 2),
      "Clover_EC", // its name
      logicEndCap, // its logical volume
      physiMother, // its mother
      true, // no boolean operat
      copyNo, // copy number
      true); // overlap check

  physiVacuum = new G4PVPlacement(0, // rotation
                                  G4ThreeVector(0. * mm, 0. * mm, vacuum_PosZ),
                                  "Clover_Vac", // its name
                                  logicVacuum, // its logical volume
                                  physiEndCap, // its mother
                                  true, // no boolean operat
                                  copyNo, // copy number
                                  true); // overlap check

  // define these here
  G4VisAttributes *visAttAlCap = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
  visAttAlCap->SetVisibility(true);
  visAttAlCap->SetForceWireframe(true);

  G4VisAttributes *visAttGeVac = new G4VisAttributes(G4Colour(0.9, 1.0, 0.9));
  visAttGeVac->SetForceWireframe(true);
  visAttGeVac->SetVisibility(true);

  G4VisAttributes *visAttActive = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0));
  // visAttActive->SetForceWireframe(true);
  visAttActive->SetVisibility(true);

  G4VisAttributes *visAttPassive = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0));
  visAttPassive->SetForceWireframe(true);
  visAttPassive->SetVisibility(true);

  G4VisAttributes *visAttLiContact =
      new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  visAttLiContact->SetVisibility(true);

  G4VisAttributes *visAttHole = new G4VisAttributes(G4Colour(1.0, 0.0, 1.0));
  visAttHole->SetVisibility(true);

  logicEndCap->SetVisAttributes(visAttAlCap);
  logicVacuum->SetVisAttributes(visAttGeVac);

  // Now for the placement of the leaves in each clover......
  G4RotationMatrix *rmC;
  G4double dPos = fGeLeaf_dX + fGapBetweenLeaves / 2.;
  G4double leafX;
  G4double leafY;
  G4double leafZ;

  for (G4int l = 0; l < 4; l++) {
    // the rotation
    rmC = new G4RotationMatrix;
    rmC->set(0, 0, 0);
    rmC->rotateZ(90. * degree * (4 - l));
    rmC->invert();
    // the x-translation
    if (l < 2) {
      leafX = dPos;
    } else {
      leafX = -dPos;
    }
    // the y-translation
    if (l == 0 || l == 3) {
      leafY = dPos;
    } else {
      leafY = -dPos;
    }
    // the z-translation
    leafZ = geLeaf_PosZ;

    physiGeLeaf[l] = new G4PVPlacement(rmC, // rotation
                                       G4ThreeVector(leafX, leafY, leafZ),
                                       "Clover", // its name
                                       logicGeLeaf[l], // its logical volume
                                       physiVacuum, // its mother
                                       true, // no boolean operat
                                       copyNo * 4 + l, // copy number
                                       true); // overlap check

    physiPassivated[l] = new G4PVPlacement(
        0, // rotation
        G4ThreeVector(-fHole_dX, -fHole_dY, fContact_dZ), "GePassivated",
        logicPassivated[l], physiGeLeaf[l], false, copyNo * 4 + l, true);

    physiContact[l] = new G4PVPlacement(
        0, // rotation
        G4ThreeVector(0. * mm, 0. * mm, 0.0 * mm), //-fContact_dZ),
        "LiContact", logicContact[l], physiPassivated[l], false, copyNo * 4 + l,
        true);

    physiBoreHole[l] = new G4PVPlacement(
        0, // rotation
        G4ThreeVector(0. * mm, 0. * mm, 0.0 * mm), //-fContact_dZ),
        "BoreHole", logicBoreHole[l], physiContact[l], false, copyNo * 4 + l,
        true);

    // visual attributes
    logicGeLeaf[l]->SetVisAttributes(visAttActive);
    logicPassivated[l]->SetVisAttributes(visAttPassive);
    logicContact[l]->SetVisAttributes(visAttLiContact);
    logicBoreHole[l]->SetVisAttributes(visAttHole);
  }
}
