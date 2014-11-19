#include <iostream>
#include <fstream>

using namespace std;

int main(){

 char inFileName[15];
 cout << "Filename to convert: ";
 cin >> inFileName;
 ifstream inFile;
 inFile.open(inFileName);
 if (!inFile) 
  { cerr << "Can't open input file " << endl;
    return 1;
  }

 char outFileName[15];
 cout << "Output filename: ";
 cin >> outFileName;
 ofstream outFile;
 outFile.open(outFileName, ios::out);

 double intensity =1.0;
 int e1,e2,e3,e4,e5,e6,e7,e8,e9,e10;

 while( intensity > 0.0 )
  {
   inFile >> intensity >> e1 >> e2 >> e3 >> e4 >> e5 
                       >> e6 >> e7 >> e8 >> e9 >> e10;
   
   if (intensity > 0.00001)
    {
     outFile << intensity*100.0 << " " << e1 << " " << e2 << " " << e3
                                << " " << e4 << " " << e5 << " " << e6 
                                << " " << e7 << " " << e8 << " " << e9 
                                << " " << e10 << endl;
    }

  }

 inFile.close();
 outFile.close();

 return 0;

}


