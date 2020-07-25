#include "Structure.h"

Structure structure;
string infilename; 
string outfilename;
    
int main(int argc, char *argv[]){

  cout << endl ;
  cout << "Welcome to xtl2cif, which makes a cif from an xtl file." << endl;
  cout << "Written by Stephen Wells." << endl;
  cout << endl;
  
  if ( argc < 3  ){
    cout << "Usage:" << endl;
    cout << argv[0] << " NameIn NameOut" << endl;
    cout << "Reads structure from XTL file NameIn and output CIF to NameOut" << endl;
    cout << "Conversion is utterly mindless and reports P1 symmetry." << endl;
    exit(0) ;
  }
  else {
    
  cout << "Will read from file: " << argv[1] << endl;
  infilename = argv[1];     

  cout << "Will write to file: " << argv[2] << endl;
  outfilename = argv[2];
  }


  if ( infilename == outfilename ){
    cout << "Don't use same name for input and output!" << endl;
    exit(1);     
  }


  structure.readxtl( infilename );
  // input atoms now in atom vector
  cerr << "I have read " << structure.atom.size()  << " atoms." << endl;
  structure.writecif( outfilename, infilename ); // records infilename in title

  return(0);
} //ENDOFMAIN
