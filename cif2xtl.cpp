#include "Structure.h"

Structure structure;
string infilename; 
string outfilename;
    
int main(int argc, char *argv[]){

  cout << endl ;
  cout << "Welcome to cif2xtl, which makes an xtl from a cif file." << endl;
  cout << "Written by Stephen Wells." << endl;
  cout << endl;
  
  if ( argc < 3  ){
    cout << "Usage:" << endl;
    cout << argv[0] << " NameIn NameOut" << endl;
    cout << "Reads structure from CIF file NameIn and outputs XTL to NameOut" << endl;
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


  structure.readcif( infilename );
  // input atoms now in atom vector
  cerr << "I have read " << structure.atom.size()  << " atoms." << endl;
  structure.writextl( outfilename, infilename ); // records infilename in title

  return(0);
} //ENDOFMAIN
