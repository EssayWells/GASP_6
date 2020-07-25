#include "Structure.h"

vector< string > header;

Structure structure;

int nx, ny, nz;

string infilename; 
string outfilename;

int main(int argc, char *argv[]){

  cout << endl << endl; 
  cout << "Welcome to Super, which makes a supercell from a CIF file." << endl;
  cout << "Written by Stephen Wells." << endl;
  cout << endl;

  if ( argc < 6  ){
    cout << "Usage:" << endl;
    cout << argv[0] << " NameIn NX NY NZ NameOut" << endl;
    cout << "Reads structure from CIF file NameIn and output to NameOut" << endl;
    cout << "Creates a supercell by multiplying up by NX, NY, NZ (integers!)" << endl;
    exit(0) ;
  }
  else {
    
  cout << "Will read from file: " << argv[1] << endl;
  infilename = argv[1];     

  cout << "Will write to file: " << argv[5] << endl;
  outfilename = argv[5];
  }



  if ( infilename == outfilename ){
    cout << "Don't use same name for input and output!" << endl;
    exit(1);     
  }

  nx = atoi( argv[2] );
  ny = atoi( argv[3] );
  nz = atoi( argv[4] );  

  cout << "Supercell size: " << nx << "," << ny << "," << nz << endl;


  structure.readcif( infilename );
  // input atoms now in atom vector
  structure.super( nx, ny, nz );

  stringstream otitle;
  otitle << " supercell from " << infilename << " " << nx << " " << ny << " " << nz;
  structure.writecif( outfilename, otitle.str() );

  return(0);
} //ENDOFMAIN
