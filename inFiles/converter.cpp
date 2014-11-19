
#include <iostream>
#include <fstream>

using namespace std;

int main(){

 //Get input file
   char inFileName[15];
   cout << "Filename to convert: ";
   cin >> inFileName;
   ifstream inFile;
   inFile.open(inFileName);
   if (!inFile) 
    { cerr << "Can't open input file " << endl;
      return 1;
    }

 //Initialize output file
   char outFileName[15];
   cout << "Output filename: ";
   cin >> outFileName;
   ofstream outFile;
   outFile.open(outFileName, ios::out);

 //Variables
   double intensity =1.0; //gamma cascade intensity from input file
   int eIn[10];           //gamma cascade energies from input file
   int eOut[10];          //energies to output to geant4 file
   int multiplicity;      //multiplicity of cascade


 //Initialize output energies to 0
   for (int i=0; i<10; i++)
    {eOut[i]=0;}

 
 while( intensity > 0.0 )
  {
   inFile >> intensity;    //read in intensity 
   for(int i=0; i<10; i++) //read in energies
    {inFile >> eIn[i];}    

   //calculate multiplicity
     if ((eIn[0]>0) && (eIn[1]==0)) {multiplicity =1;}
     if ((eIn[1]>0) && (eIn[2]==0)) {multiplicity =2;}
     if ((eIn[2]>0) && (eIn[3]==0)) {multiplicity =3;}
     if ((eIn[3]>0) && (eIn[4]==0)) {multiplicity =4;}
     if ((eIn[4]>0) && (eIn[5]==0)) {multiplicity =5;}
     if ((eIn[5]>0) && (eIn[6]==0)) {multiplicity =6;}
     if ((eIn[6]>0) && (eIn[7]==0)) {multiplicity =7;}
     if ((eIn[7]>0) && (eIn[8]==0)) {multiplicity =8;}
     if ((eIn[8]>0) && (eIn[9]==0)) {multiplicity =9;}
     if (eIn[9]>0) {multiplicity =10;}

   //Only keep intensities that are big enough   
     if (intensity > 0.00001)
      {
       //write out intensity
       outFile << intensity*100.0 << "       " ;

       //figure out the output energy that we want     
       for(int i=0; i<multiplicity; i++)
        {
          if(eIn[i]>0)
          eOut[i]=eIn[i];
        }
       for(int i=multiplicity; i<10; i++)
        {
          eOut[i]=0;      
        }

  
     //write out energies 
       for(int i=0; i<10; i++)
        {
          outFile << eOut[i] << " " ;
        }
        outFile << endl;
      }

  }


 //Close files
   inFile.close();
   outFile.close();

 return 0;

}


