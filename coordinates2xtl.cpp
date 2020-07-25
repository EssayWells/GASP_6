#include "Structure.h" // uses GASP Structure header, which brings in Rotor, Vector, basic_utils

Structure structure;
string infilename; 
string outfilename;

void read_coords( string filename);

int main(int argc, char *argv[]){

  cerr << endl << endl; 
  cerr << "Welcome to coordinates2xtl, which makes an xtl file from a Crystalmaker coordinates.txt file." << endl;
  cerr << "Written by Stephen Wells." << endl;
  cerr << endl;

  if ( argc < 3  ){
    cerr << "Usage:" << endl;
    cerr << argv[0] << " NameIn NameOut" << endl;
    cerr << "Reads structure from Coordinates file NameIn and outputs xtl to NameOut" << endl;
    cerr << "Conversion is utterly mindless and reports P1 symmetry." << endl;
    exit(0) ;
  }
  else {
    
	cerr << "Will read from file: " << argv[1] << endl;
	infilename = argv[1];     

	cerr << "Will write to file: " << argv[2] << endl;
	outfilename = argv[2];
  }

  if ( infilename == outfilename ){
    cerr << "Don't use same name for input and output!" << endl;
    exit(1);     
  }


  bool goodin = structure.read_coords( infilename );
  if ( goodin ){
	  cerr << "Read coordinates ok." << endl;
  }
  else{
	  cerr << "Bad coordinates file?" << endl;
  }
  // input atoms now in structure

  stringstream outtitle;
  outtitle << ": xtl from " << infilename;
  structure.writextl( outfilename, outtitle.str() , false );

  return(0);
} //ENDOFMAIN

