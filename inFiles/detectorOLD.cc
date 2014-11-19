// SuN = Summing NaI detector
 
#include "XriDetectorConstruction.hh"
#include "G4SDManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4UnitsTable.hh"
#include "globals.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <iostream>
#include <sstream>
#include "G4String.hh"
#include "G4ios.hh"
#include <stdio.h>

XriDetectorConstruction::XriDetectorConstruction()
  :  NaI(0), Al(0), N78O21Ar1(0), Cr20Ni8Fe76(0), O65Si25B7Na2Al1(0)
{
}

XriDetectorConstruction::~XriDetectorConstruction()
{
}

G4VPhysicalVolume* XriDetectorConstruction::Construct()
{
  DefineMaterials();

 return ConstructDetector();
}

void XriDetectorConstruction::DefineMaterials()
{
 //This function illustrates the possible ways to define materials

	G4String name, symbol;             //a=mass of a mole;
	G4double a, z, density;            //z=mean number of protons;

	G4int ncomponents, natoms;
      //G4double abundance, fractionmass;
      //G4double temperature, pressure;



// define Elements

	a = 22.99*g/mole;
	G4Element* Na = new G4Element(name="Sodium" ,symbol="Na" , z= 11., a);

	a = 126.90*g/mole;
	G4Element* I = new G4Element(name="Iodine" ,symbol="I" , z= 53., a);

	a = 204.38*g/mole;
	G4Element* Tl = new G4Element(name="Thalium" ,symbol="Tl" , z= 81., a);

	a = 26.982*g/mole;
	G4Element* elAl  = new G4Element(name="element_Aluminum",symbol="elAl" , z= 13., a);

	a = 14.00*g/mole;
	G4Element* N  = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);

	a = 16.00*g/mole;
	G4Element* O  = new G4Element(name="Oxygen",symbol="O" , z= 8., a);

	a = 39.95*g/mole;
	G4Element* Ar  = new G4Element(name="Argon",symbol="Ar" , z= 18., a);

        a = 51.996*g/mole;
	G4Element* Cr  = new G4Element(name="Chromium"  ,symbol="Cr" , z= 24., a);

        a = 58.69*g/mole;
	G4Element* Ni  = new G4Element(name="Nickel" ,symbol="Ni" , z= 28., a);

	a = 55.847*g/mole;
	G4Element* Fe  = new G4Element(name="Iron"  ,symbol="Fe" , z= 26., a);
	    
	a = 28.086*g/mole;
	G4Element* Si = new G4Element(name="Silicon" ,symbol="Si" , z= 14., a);

	a = 10.811*g/mole;
	G4Element* B = new G4Element(name="Boron" ,symbol="B" , z= 5., a);


//____________________________________________________________________________
// define a material from elements.   case 1: chemical molecule


// ........Stainless Steel..................

	density = 8.0*g/cm3;
	Cr20Ni8Fe76 = new G4Material(name="Stainless_Steel", density, ncomponents=3);
	Cr20Ni8Fe76->AddElement(Cr, natoms=20);
        Cr20Ni8Fe76->AddElement(Fe, natoms=76);
	Cr20Ni8Fe76->AddElement(Ni, natoms=8);

// ...........Pyrex........................

	density = 2.23*g/cm3;
	O65Si25B7Na2Al1 = new G4Material(name="Pyrex", density, ncomponents=5);
	O65Si25B7Na2Al1->AddElement(O, natoms=65);
        O65Si25B7Na2Al1->AddElement(Si, natoms=25);
	O65Si25B7Na2Al1->AddElement(B, natoms=7);
        O65Si25B7Na2Al1->AddElement(Na, natoms=2);
        O65Si25B7Na2Al1->AddElement(elAl, natoms=1);

//.............NaI...........................
	
	density = 3.67*g/cm3;
	NaI = new G4Material(name="Sodium Iodine", density, ncomponents=3);
	NaI->AddElement(Na, natoms=1000);
	NaI->AddElement(I, natoms=1000);
        NaI->AddElement(Tl, natoms=1);

//................Al...............	
	
	density = 2.698*g/cm3;
	Al = new G4Material (name="Aluminum", density, ncomponents=1);
	Al->AddElement(elAl, natoms=1);
	
	
//.............Air......................

	density = 1.2927*mg/cm3;
	N78O21Ar1 = new G4Material (name="Air", density, ncomponents=3);
	N78O21Ar1->AddElement(N, natoms=78);
	N78O21Ar1->AddElement(O, natoms=21);
	N78O21Ar1->AddElement(Ar, natoms=1);
	
	
//........................................................................
  G4cout << "\n\n ####-------------------------------------------------------#### \n";
  G4cout << "\n\n\n\n\t\t #### List of elements used #### \n";
  G4cout << *(G4Element::GetElementTable());
  G4cout << "\n\n\n\n\t\t #### List of materials used #### \n";
  G4cout << *(G4Material::GetMaterialTable());
  G4cout << "\n\n ####-------------------------------------------------------#### \n";

}

G4VPhysicalVolume* XriDetectorConstruction::ConstructDetector()
{


//..........EXPERIMENTAL ROOM........
  G4Tubs* room_tube = new G4Tubs("room", 0.0*cm, 100.0*cm, 300.0*cm, 0.0*deg, 360.0*deg);
  G4LogicalVolume* room_log = new G4LogicalVolume(room_tube,N78O21Ar1,"room",0,0,0);
  G4VPhysicalVolume* room_phys = new G4PVPlacement(0,G4ThreeVector(0.0*cm,0.0*cm,0.0*cm),"room",room_log,NULL,false,0);


//..........BEAM PIPE................
  G4double outerR_beam = 21.50*mm;                  //edit this to change the radius of the beam pipe
  G4double innerR_beam = outerR_beam - 1.0*mm;      //edit this to change the thickness of the beam pipe
  G4double halflength_beam = 60.0*cm;               //edit this to change the length of the beam pipe
  G4double startAngle_beam = 0.*deg;
  G4double spanAngle_beam = 360.*deg;

  G4Tubs* beam_tube = new G4Tubs("beam_tube",innerR_beam, outerR_beam, halflength_beam, startAngle_beam, spanAngle_beam);

  G4LogicalVolume* beam_log = new G4LogicalVolume(beam_tube,Cr20Ni8Fe76,"beam_log",0,0,0);

  G4VPhysicalVolume* beam_phys = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,0.0*mm),beam_log,"beam_phys",room_log,false,0);


//........NaI SCINTILLATOR........
  G4double innerR_scint = 45.0*mm*0.5;              // 45mm borehole
  G4double outerR_scint = 203.0*mm;                 // total of 406mm in diameter
  G4double length_scint = 101.5*mm;                 // total of 406mm in length
  G4double halflength_scint = 0.5*length_scint;
  G4double startAngle_scint = 0.0*deg;
  G4double spanAngle_scint = 180.0*deg;

  G4Tubs* scint_tube = new G4Tubs("scint_tube",innerR_scint,outerR_scint,halflength_scint,startAngle_scint,spanAngle_scint);
 
  G4LogicalVolume* scint_log = new G4LogicalVolume(scint_tube,NaI,"scint_log",0,0,0);


//..........REFLECTOR............
  G4double width_Refl = 0.25*mm;               
  G4double halfwidth_Refl = 0.5*width_Refl; 

 // 1 = outside of scintillator
  G4double innerR_Refl1 = outerR_scint;
  G4double outerR_Refl1 = innerR_Refl1 + width_Refl;

  G4Tubs* reflector_tube1 = new G4Tubs("reflector_tube1",innerR_Refl1,outerR_Refl1,halflength_scint,startAngle_scint,spanAngle_scint);
 
  G4LogicalVolume*  reflector_log1 = new G4LogicalVolume(reflector_tube1,O65Si25B7Na2Al1,"reflector_log1",0,0,0);

 // 2 = inside of scintillator
   G4double outerR_Refl2 = innerR_scint;
   G4double innerR_Refl2 = outerR_Refl2 - width_Refl;

   G4Tubs* reflector_tube2 = new G4Tubs("reflector_tube2",innerR_Refl2,outerR_Refl2,halflength_scint,startAngle_scint,spanAngle_scint);
  
   G4LogicalVolume*  reflector_log2 = new G4LogicalVolume(reflector_tube2,O65Si25B7Na2Al1,"reflector_log2",0,0,0);

 // 3 = side of scintillator 
   G4double innerR_Refl3 = innerR_Refl2;
   G4double outerR_Refl3 = outerR_Refl1;

   G4Tubs* reflector_tube3 = new G4Tubs("reflector_tube3",innerR_Refl3,outerR_Refl3,halfwidth_Refl,startAngle_scint,spanAngle_scint);
 
   G4LogicalVolume*  reflector_log3 = new G4LogicalVolume(reflector_tube3,O65Si25B7Na2Al1,"reflector_log3",0,0,0);

 // 4 = flat middle sections 
  G4double halfXlength_Refl4 = (outerR_Refl1-innerR_Refl2)*0.5;
  G4double halfYlength_Refl4 = halfwidth_Refl;
  G4double halfZlength_Refl4 = halflength_scint + width_Refl;

  G4Box *reflector_box = new G4Box("reflector_box",halfXlength_Refl4,halfYlength_Refl4,halfZlength_Refl4);

  G4LogicalVolume *reflector_log4 = new G4LogicalVolume(reflector_box,O65Si25B7Na2Al1,"reflector_log4",0,0,0);



//..........ALUMINUM.............
  G4double width_Al = 0.50*mm;
  G4double halfwidth_Al = 0.5*width_Al;

 // 1 = outside of scintillator
  G4double innerR_Al1 = outerR_Refl1;
  G4double outerR_Al1 = innerR_Al1 + width_Al;
  G4double halflength_Al1 = halflength_scint + 2*width_Refl;

  G4Tubs* Al_tube1 = new G4Tubs("Al_tube1",innerR_Al1,outerR_Al1,halflength_Al1,startAngle_scint,spanAngle_scint);
 
  G4LogicalVolume*  Al_log1 = new G4LogicalVolume(Al_tube1,Al,"Al_log1",0,0,0);

 // 2 = inside of scintillator
  G4double outerR_Al2 = innerR_Refl2;
  G4double innerR_Al2 = outerR_Al2 - width_Al;
  G4double halflength_Al2 = halflength_scint + 2*width_Refl;

  G4Tubs* Al_tube2 = new G4Tubs("Al_tube2",innerR_Al2,outerR_Al2,halflength_Al2,startAngle_scint,spanAngle_scint);
 
  G4LogicalVolume*  Al_log2 = new G4LogicalVolume(Al_tube2,Al,"Al_log2",0,0,0);

 // 3 = side of scintillator 
  G4double innerR_Al3 = innerR_Al2;
  G4double outerR_Al3 = outerR_Al1;

  G4Tubs* Al_tube3 = new G4Tubs("Al_tube3",innerR_Al3,outerR_Al3,halfwidth_Al,startAngle_scint,spanAngle_scint);
 
  G4LogicalVolume*  Al_log3 = new G4LogicalVolume(Al_tube3,Al,"Al_log3",0,0,0);

 // 4 = flat middle sections 
  G4double halfXlength_Al4 = (outerR_Al1-innerR_Al2)*0.5;
  G4double halfYlength_Al4 = halfwidth_Al;
  G4double halfZlength_Al4 = halflength_scint+width_Refl+width_Al;

  G4Box *Al_box1 = new G4Box("Al_box1",halfXlength_Al4,halfYlength_Al4,halfZlength_Al4);

  G4LogicalVolume *Al_log4 = new G4LogicalVolume(Al_box1,Al,"Al_log4",0,0,0);

 // 5,6 = filling in the gaps    
  G4double halfXlength_Al5 = 0.5*(outerR_Al1-innerR_Al2);
  G4double halfYlength_Al5 = halfwidth_Refl;
  G4double halfZlength_Al5 = halfwidth_Al;

  G4Box *Al_box2 = new G4Box("Al_box2",halfXlength_Al5,halfYlength_Al5,halfZlength_Al5);

  G4LogicalVolume *Al_log5 = new G4LogicalVolume(Al_box2,Al,"Al_log5",0,0,0);

  G4double halfXlength_Al6 = halfwidth_Al;
  G4double halfYlength_Al6 = halfwidth_Refl;
  G4double halfZlength_Al6 = halflength_scint+width_Refl+width_Al;

  G4Box *Al_box3 = new G4Box("Al_box3",halfXlength_Al6,halfYlength_Al6,halfZlength_Al6);

  G4LogicalVolume *Al_log6 = new G4LogicalVolume(Al_box3,Al,"Al_log6",0,0,0);

 // 7 = thick outer casing
  G4double innerR_Al7 = 0.0*mm;
  G4double outerR_Al7 = outerR_Al1;
  G4double width_Al7 = 12.0*mm;
  G4double halfwidth_Al7 = 0.5 * width_Al7;

  G4Tubs* Al_tube7 = new G4Tubs("Al_tube7",innerR_Al7,outerR_Al7,halfwidth_Al7,0.0*deg,360.0*deg);
 
  G4LogicalVolume*  Al_log7 = new G4LogicalVolume(Al_tube7,Al,"Al_log7",0,0,0);


  



//_______________the rotation matrix________________
  G4RotationMatrix* rot_180 = new G4RotationMatrix();
    rot_180->rotateZ(180*deg);

 

//...............NEW TARGET......................

  G4double innerR_tube1 = 3.175*0.5*mm;
  G4double outerR_tube1 = 3.81*0.5*mm;               
  G4double halflength_tube1 =15.0*cm;            
  G4double startAngle_tube1 = 0.0*deg;
  G4double spanAngle_tube1 = 360.0*deg;

  G4double innerR_tube2 = 6.35*0.5*mm;
  G4double outerR_tube2 = 6.985*0.5*mm;               
  G4double halflength_tube2 =15.25*cm;            
  G4double startAngle_tube2 = 0.0*deg;
  G4double spanAngle_tube2 = 360.0*deg;

  G4Tubs* target_tube1 = new G4Tubs("target_tube1",innerR_tube1, outerR_tube1, halflength_tube1, startAngle_tube1, spanAngle_tube1);
  G4Tubs* target_tube2 = new G4Tubs("target_tube2",innerR_tube2, outerR_tube2, halflength_tube2, startAngle_tube2, spanAngle_tube2);

  G4LogicalVolume* target_log1 = new G4LogicalVolume(target_tube1,Al,"target_log1",0,0,0);
  G4LogicalVolume* target_log2 = new G4LogicalVolume(target_tube2,Al,"target_log2",0,0,0);

  G4VPhysicalVolume* target_phys1r = new G4PVPlacement(0,G4ThreeVector(17.0*mm,0.*mm,halflength_tube1+0.5*cm),target_log1,"target_phys1r",room_log,false,0);

  G4VPhysicalVolume* target_phys2r = new G4PVPlacement(0,G4ThreeVector(17.0*mm,0.*mm,halflength_tube2),target_log2,"target_phys2r",room_log,false,0);

  G4VPhysicalVolume* target_phys1l = new G4PVPlacement(0,G4ThreeVector(-17.0*mm,0.*mm,halflength_tube1+0.5*cm),target_log1,"target_phys1l",room_log,false,0);

  G4VPhysicalVolume* target_phys2l = new G4PVPlacement(0,G4ThreeVector(-17.0*mm,0.*mm,halflength_tube2),target_log2,"target_phys2l",room_log,false,0);


  G4double innerR_backing = 0.0*mm;
  G4double outerR_backing = 20.5*mm;               
  G4double halflength_backing =2.7*mm;            
  G4double startAngle_backing = 0.0*deg;
  G4double spanAngle_backing = 360.0*deg;

  G4Tubs* backing_tube = new G4Tubs("backing_tube",innerR_backing, outerR_backing, halflength_backing, startAngle_backing, spanAngle_backing);

  G4LogicalVolume* backing_log = new G4LogicalVolume(backing_tube,Al,"backing_log",0,0,0);

  G4VPhysicalVolume* backing_phys = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,-halflength_backing),backing_log,"backing_phys",room_log,false,0);


//////////////////////////////////////////////////////
//        Placing the volumes                       //
// //////////////////////////////////////////////////


  G4double Pos_x = 0.0*cm;
  G4double Pos_y = width_Al + width_Refl;
  G4double Pos_z = -3*width_Al - 3*width_Refl -3*halflength_scint;  
  G4int inc = 0;

  detno=0;
  std::ostringstream oss;
  G4String s;

while(inc<4)
 {

 //______________ THE TOP OF THE DETECTOR ___________

   Pos_y = width_Al + width_Refl;

     ++detno;
     oss.str("");
     oss << detno;
     s=oss.str();
     detname[detno]="NaI" + s;

      X[detno]=Pos_x;
      Y[detno]=Pos_y;
      Z[detno]=Pos_z;


   //SCINTILLATOR
     G4VPhysicalVolume* scint_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z),scint_log,detname[detno],room_log,false,0);

   //REFLECTOR  
     G4VPhysicalVolume* reflector_outer_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z),reflector_log1,"reflector_outer_top",room_log,false,0);

     G4VPhysicalVolume* reflector_inner_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z),reflector_log2,"reflector_inner_top",room_log,false,0);

     G4VPhysicalVolume* reflector_sideA_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z-halflength_scint-halfwidth_Refl),reflector_log3,"reflector_sideA_top",room_log,false,0);

     G4VPhysicalVolume* reflector_sideB_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z+halflength_scint+halfwidth_Refl),reflector_log3,"reflector_sideB_top",room_log,false,0);

     G4VPhysicalVolume* reflector_middleA_top = new G4PVPlacement(0,G4ThreeVector((innerR_Refl2+outerR_Refl1)*0.5,width_Al+halfwidth_Refl,Pos_z),reflector_log4,"reflector_middleA_top",room_log,false,0);

     G4VPhysicalVolume* reflector_middleB_top = new G4PVPlacement(0,G4ThreeVector((innerR_Refl2+outerR_Refl1)*-0.5,width_Al+halfwidth_Refl,Pos_z),reflector_log4,"reflector_middleB_top",room_log,false,0);

   //ALUMINUM
     G4VPhysicalVolume* Al_outer_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z),Al_log1,"Al_outer_top",room_log,false,0);

     G4VPhysicalVolume* Al_inner_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z),Al_log2,"Al_inner_top",room_log,false,0);

     G4VPhysicalVolume* Al_sideA_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z-halflength_scint-width_Refl-halfwidth_Al),Al_log3,"Al_sideA_top",room_log,false,0);

     G4VPhysicalVolume* Al_sideB_top = new G4PVPlacement(0,G4ThreeVector(Pos_x,Pos_y,Pos_z+halflength_scint+width_Refl+halfwidth_Al),Al_log3,"Al_sideB_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleA_top = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al2),halfwidth_Al,Pos_z),Al_log4,"Al_middleA_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleB_top = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al2),halfwidth_Al,Pos_z),Al_log4,"Al_middleB_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleC_top = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al2),width_Al+halfwidth_Refl,Pos_z-halflength_scint-width_Refl-halfwidth_Al),Al_log5,"Al_middleC_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleD_top = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al2),width_Al+halfwidth_Refl,Pos_z+halflength_scint+width_Refl+halfwidth_Al),Al_log5,"Al_middleD_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleE_top = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al2),width_Al+halfwidth_Refl,Pos_z-halflength_scint-width_Refl-halfwidth_Al),Al_log5,"Al_middleE_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleF_top = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al2),width_Al+halfwidth_Refl,Pos_z+halflength_scint+width_Refl+halfwidth_Al),Al_log5,"Al_middleF_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleG_top = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al1),width_Al+halfwidth_Refl,Pos_z),Al_log6,"Al_middleG_top",room_log,false,0);

     G4VPhysicalVolume* Al_middleH_top = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al1),width_Al+halfwidth_Refl,Pos_z),Al_log6,"Al_middleH_top",room_log,false,0);




 //____________ THE BOTTOM OF THE DETECTOR ___________________
 
   Pos_y = -width_Al - width_Refl;

     ++detno;
     oss.str("");
     oss << detno;
     s=oss.str();
     detname[detno]="NaI" + s;

      X[detno]=Pos_x;
      Y[detno]=Pos_y;
      Z[detno]=Pos_z;

   //SCINTILLATOR
     G4VPhysicalVolume* scint_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z),scint_log,detname[detno],room_log,false,0);

   //REFLECTOR  
     G4VPhysicalVolume* reflector_outer_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z),reflector_log1,"reflector_outer_bottom",room_log,false,0);

     G4VPhysicalVolume* reflector_inner_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z),reflector_log2,"reflector_inner_bottom",room_log,false,0);

     G4VPhysicalVolume* reflector_sideA_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z-halflength_scint-halfwidth_Refl),reflector_log3,"reflector_sideA_bottom",room_log,false,0);

     G4VPhysicalVolume* reflector_sideB_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z+halflength_scint+halfwidth_Refl),reflector_log3,"reflector_sideB_bottom",room_log,false,0);

     G4VPhysicalVolume* reflector_middleA_bottom = new G4PVPlacement(0,G4ThreeVector((innerR_Refl2+outerR_Refl1)*0.5,-width_Al-halfwidth_Refl,Pos_z),reflector_log4,"reflector_middleA_bottom",room_log,false,0);

     G4VPhysicalVolume* reflector_middleB_bottom = new G4PVPlacement(0,G4ThreeVector((innerR_Refl2+outerR_Refl1)*-0.5,-width_Al-halfwidth_Refl,Pos_z),reflector_log4,"reflector_middleB_bottom",room_log,false,0);

   //ALUMINUM
     G4VPhysicalVolume* Al_outer_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z),Al_log1,"Al_outer_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_inner_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z),Al_log2,"Al_inner_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_sideA_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z-halflength_scint-width_Refl-halfwidth_Al),Al_log3,"Al_sideA_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_sideB_bottom = new G4PVPlacement(rot_180,G4ThreeVector(Pos_x,Pos_y,Pos_z+halflength_scint+width_Refl+halfwidth_Al),Al_log3,"Al_sideB_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleA_bottom = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al2),-halfwidth_Al,Pos_z),Al_log4,"Al_middleA_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleB_bottom = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al2),-halfwidth_Al,Pos_z),Al_log4,"Al_middleB_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleC_bottom = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al2),-width_Al-halfwidth_Refl,Pos_z-halflength_scint-width_Refl-halfwidth_Al),Al_log5,"Al_middleC_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleD_bottom = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al2),-width_Al-halfwidth_Refl,Pos_z+halflength_scint+width_Refl+halfwidth_Al),Al_log5,"Al_middleD_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleE_bottom = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al2),-width_Al-halfwidth_Refl,Pos_z-halflength_scint-width_Refl-halfwidth_Al),Al_log5,"Al_middleE_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleF_bottom = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al2),-width_Al-halfwidth_Refl,Pos_z+halflength_scint+width_Refl+halfwidth_Al),Al_log5,"Al_middleF_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleG_bottom = new G4PVPlacement(0,G4ThreeVector(0.5*(outerR_Al1+innerR_Al1),-width_Al-halfwidth_Refl,Pos_z),Al_log6,"Al_middleG_bottom",room_log,false,0);

     G4VPhysicalVolume* Al_middleH_bottom = new G4PVPlacement(0,G4ThreeVector(-0.5*(outerR_Al1+innerR_Al1),-width_Al-halfwidth_Refl,Pos_z),Al_log6,"Al_middleH_bottom",room_log,false,0);

   Pos_z = Pos_z + 2*width_Al + 2*width_Refl + length_scint;
   ++inc;

   }

  Pos_z = 2.0*length_scint + 4.0*width_Al + 4.0*width_Refl + halfwidth_Al7; 
  G4VPhysicalVolume* Al_thickA = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,Pos_z),Al_log7,"Al_thickA",room_log,false,0);
  G4VPhysicalVolume* Al_thickB = new G4PVPlacement(0,G4ThreeVector(0.0*mm,0.0*mm,-Pos_z),Al_log7,"Al_thickB",room_log,false,0);




//========================== Visualization attributes =========================================//

  room_log->SetVisAttributes (G4VisAttributes::Invisible);

//visualization for scintillators =red
  G4VisAttributes *RedAttr = new G4VisAttributes(G4Colour(1.,0.,0.));      //Red
  RedAttr->SetVisibility(true);
  // RedAttr->SetForceWireframe(true);
  RedAttr->SetForceSolid(true);

//visualization for aluminum = green
  G4VisAttributes *GreenAttr = new G4VisAttributes(G4Colour(0.,1.,0.));      //Green
  GreenAttr->SetVisibility(true);
  // GreenAttr->SetForceWireframe(true);
  GreenAttr->SetForceSolid(true);

//visualization for beam = grey
  G4VisAttributes *GreyAttr = new G4VisAttributes(G4Colour(0.5,0.5,0.5));      //Grey
  GreyAttr->SetVisibility(true);
  // ScintubeAttr->SetForceWireframe(true);
  GreyAttr->SetForceSolid(true);

//visualization for blue
  G4VisAttributes *BlueAttr = new G4VisAttributes(G4Colour(0.,0.,1.));      //Blue
  BlueAttr->SetVisibility(true);
  // BlueAttr->SetForceWireframe(true);
  BlueAttr->SetForceSolid(true);

//visualization for purple
  G4VisAttributes *PurpleAttr = new G4VisAttributes(G4Colour(1.,0.,1.));     //Purple
  PurpleAttr->SetVisibility(true);
  // PurpleAttr->SetForceWireframe(true);
  PurpleAttr->SetForceSolid(true);


 //applying the color scheme

scint_log->SetVisAttributes(RedAttr);
beam_log->SetVisAttributes(GreyAttr);
reflector_log1->SetVisAttributes(BlueAttr);
reflector_log2->SetVisAttributes(BlueAttr);
reflector_log3->SetVisAttributes(BlueAttr);
reflector_log4->SetVisAttributes(BlueAttr);
Al_log1->SetVisAttributes(GreenAttr);
Al_log2->SetVisAttributes(GreenAttr);
Al_log3->SetVisAttributes(GreenAttr);
Al_log4->SetVisAttributes(GreenAttr);
Al_log5->SetVisAttributes(GreenAttr);  
Al_log6->SetVisAttributes(GreenAttr);
Al_log7->SetVisAttributes(GreenAttr);

target_log1->SetVisAttributes(GreyAttr);
target_log2->SetVisAttributes(GreyAttr);
backing_log->SetVisAttributes(GreyAttr);

  return room_phys;
}



