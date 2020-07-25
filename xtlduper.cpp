#include "Structure.h"

vector< string > header;

Structure structure; // read in; then strip
//so give Structure a new strip-dupes power!

int main(int argc, char *argv[]){

  cout << endl << endl; 
  cout << "Welcome to Duper, which cuts duplicate atoms from xtl files." << endl;
  cout << "Written by Stephen Wells." << endl;
  cout << endl;
  
    string infilename; 
    string outfilename;
  
  if ( argc != 3  ){
    cout << "Usage:" << endl;
    cout << argv[0] << " name1 name2 : read structure from XTL file name1 and output to name2" << endl;
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
  // input atoms now in structure.atom vector
  cout << "Found " << structure.atom.size() << " input atoms." << endl;
  
  structure.stripdupes();
  cout << "Found " << structure.atom.size() << " unduplicated atoms." << endl;

  stringstream otitle;
  otitle << " from " << outfilename;
  structure.writextl( outfilename, otitle.str(), false ); // using structure write out method

  return(0);
} //ENDOFMAIN
