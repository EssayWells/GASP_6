#include "Structure.h" // brings all utilities

double dx, dy, dz;
string infilename; 
string outfilename;
Structure structure;

int main(int argc, char *argv[]){

  cout << endl << endl; 
  cout << "Welcome to shifter, which shifts all the fractional coordinates in an xtl file." << endl;
  cout << "Written by Stephen Wells." << endl;
  cout << endl;
  
  if ( argc < 6  ){
    cout << "Usage:" << endl;
    cout << argv[0] << " NameIn dX dY dZ NameOut" << endl;
    cout << "Reads structure from XTL file NameIn and outputs to NameOut" << endl;
    cout << "Increments all fractional coordinates by dX, dY, dZ and wraps into the cell again." << endl;
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

  
  dx = atof( argv[2] );
  dy = atof( argv[3] );
  dz = atof( argv[4] );  
  Vector move;
  move = Vector( dx, dy, dz );

  cout << "Required shift: " << move << endl;

  structure.readxtl( infilename );
  structure.shift( move );
  structure.writextl( outfilename, infilename, false );

  return(0);
} //ENDOFMAIN
