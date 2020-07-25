#include "Structure.h"

Structure structure1,structure2;
string infilename1; 
string infilename2;
    
int main(int argc, char *argv[]){

  cout << endl ;
  cout << "Welcome to netdr2, computes the net difference between two structures." << endl;
  cout << "Written by Stephen Wells." << endl;
  cout << endl;
  
  if ( argc < 3  ){
    cout << "Usage:" << endl;
    cout << argv[0] << " Name1 Name2" << endl;
    cout << "Reads structures from CIF files Name1 and Name2." << endl;
    cout << "The user should ensure they have matching atoms (number and order)," << endl;
    cout << "and that the cells are at least very similar if not identical." << endl;
    cout << "Will report the mean squared change in fractional coordinates, " << endl;
    cout << "with and without a net global motion correction," << endl;
    cout << "and the equivalent mean squared motion in Angstroms-squared, " << endl;
    cout << "computed using the cell parameters of file Name1. " << endl;
    cout << "Output will go to screen through the cout channel." << endl;
    exit(0) ;
  }
  else {
    
  cout << "Will read from file: " << argv[1] << endl;
  infilename1 = argv[1];     

  cout << "Will read from file: " << argv[2] << endl;
  infilename2 = argv[2];
  }


  if ( infilename1 == infilename2 ){
    cout << "Don't use same name for noth!" << endl;
    exit(1);     
  }


  structure1.readcif( infilename1 );
  // input atoms now in atom vector
  cerr << "I have read " << structure1.atom.size()  << " atoms." << endl;
  structure2.readcif( infilename2 );
  // input atoms now in atom vector
  cerr << "I have read " << structure2.atom.size()  << " atoms." << endl;  
  
  if ( structure1.atom.size() != structure2.atom.size() ){
	  cerr << "Error, number of atoms not matched." << endl;
	  exit(1);
  }
  
  //set up cells
  structure1.cell.makecell();
  structure1.cell.makestar();
  
  Vector COMdiff;
  vector< Vector > fdiff, cdiff;
  
  for ( int i = 0; i < structure1.atom.size() ; i++ ){
	  Vector adiff = structure1.atom.at(i).fracpos - structure2.atom.at(i).fracpos ;
	  //wrap!
	  Vector bdiff = wrapfracdel( adiff );
	  fdiff.push_back( bdiff );
	  COMdiff += bdiff; // track COM value
  }
  COMdiff /= structure1.atom.size(); //average it out
  cdiff = fdiff;
  double netUncorrectedFrac = 0.;
  double netCorrectedFrac = 0.;
  for ( int i = 0; i < structure1.atom.size() ; i++ ){
	  netUncorrectedFrac += fdiff.at(i).sq() ; //self-dot added
	  cdiff.at(i) -= COMdiff;
	  netCorrectedFrac += cdiff.at(i).sq();
  }  
  double netUncorrectedCart = 0.;
  double netCorrectedCart = 0.;
  for ( int i = 0; i < structure1.atom.size() ; i++ ){
	  Vector ac = structure1.cell.fractocart( fdiff.at(i) );
	  Vector bc = structure1.cell.fractocart( cdiff.at(i) );
	  netUncorrectedCart += ac.sq();
	  netCorrectedCart += bc.sq();
  }  
  
  cout << "Mean-delta-fractional-coordinates:                      " << COMdiff << endl;
  cout << "Total-squared-delta-fractional-coordinates-uncorrected: " << netUncorrectedFrac << endl;
  cout << "Total-squared-delta-fractional-coordinates-corrected:   " << netCorrectedFrac << endl;
  cout << "Total-squared-delta-Cartesian-coordinates-uncorrected:  " << netUncorrectedCart << endl;
  cout << "Total-squared-delta-Cartesian-coordinates-corrected:    " << netCorrectedCart << endl;

  return(0);
} //ENDOFMAIN
