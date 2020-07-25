#include "Structure.h" // brings all utilities

string infilename; 
string outfilename;
Structure structure;

int x1, x2, x3, x4, x5, x6, x7, x8, x9;

int main(int argc, char *argv[]){

  cout << endl << endl; 
  cout << "Welcome to " << argv[0] << " , which transforms fractional coordinates in a cif file." << endl;
  cout << "Written by Stephen Wells." << endl;
  cout << endl;
  
  if ( argc < 12  ){
    cout << "Usage:" << endl;
    cout << argv[0] << " NameIn A B C D E F G H I NameOut" << endl;
    cout << "Reads a structure from CIF file NameIn and outputs to NameOut" << endl;
    cout << "Uses a new lattice setting defined by ABCDEFGHI." << endl;
    cout << "Each of the above must be an integer between -3 and 3." << endl;
    cout << "The new [100] vector is [A B C];" << endl;
    cout << "The new [010] vector is [D E F];" << endl;
    cout << "The new [001] vector is [G H I];" << endl;
    cout << "User is responsible for sensible choices." << endl << endl;
    exit(0) ;
  }
  else {
    
  cout << "Will read from file: " << argv[1] << endl;
  infilename = argv[1];     

  cout << "Will write to file: " << argv[2] << endl;
  outfilename = argv[11];
  }

  if ( infilename == outfilename ){
    cout << "Don't use same name for input and output!" << endl;
    exit(1);     
  }
  
  x1 = atoi( argv[2] );
  x2 = atoi( argv[3] );
  x3 = atoi( argv[4] );
  x4 = atoi( argv[5] );
  x5 = atoi( argv[6] );
  x6 = atoi( argv[7] );
  x7 = atoi( argv[8] ); 
  x8 = atoi( argv[9] );  
  x9 = atoi( argv[10] );  
  
  //test for bounds here
  int p = -3;
  int q = 3;
  if ( x1 < p || x1 > q ){
	  cerr << "Value " << x1 << " out of bounds." << endl;
	  exit(1);
  }
  if ( x2 < p || x2 > q ){
	  cerr << "Value " << x2 << " out of bounds." << endl;
	  exit(1);
  }
  if ( x3 < p || x3 > q ){
	  cerr << "Value " << x3 << " out of bounds." << endl;
	  exit(1);
  }  
  if ( x4 < p || x4 > q ){
	  cerr << "Value " << x4 << " out of bounds." << endl;
	  exit(1);
  }
  if ( x5 < p || x5 > q ){
	  cerr << "Value " << x5 << " out of bounds." << endl;
	  exit(1);
  }
  if ( x6 < p || x6 > q ){
	  cerr << "Value " << x6 << " out of bounds." << endl;
	  exit(1);
  }
  if ( x7 < p || x7 > q ){
	  cerr << "Value " << x7 << " out of bounds." << endl;
	  exit(1);
  }
  if ( x8 < p || x8 > q ){
	  cerr << "Value " << x8 << " out of bounds." << endl;
	  exit(1);
  }
  if ( x9 < p || x9 > q ){
	  cerr << "Value " << x9 << " out of bounds." << endl;
	  exit(1);
  }

  //
  
  structure.readcif( infilename );

  vector< Atom > bigatom; //build a huge cell here
  vector< Atom > smallatom; //build the single cell here
  
  int nat = structure.atom.size(); //how many atoms originally?
  cerr << "Working on " << nat << " atoms in input structure." << endl;

  //small cell first
  for ( int i = 0 ; i < nat ; i++ ){
	  Vector p = structure.atom.at(i).fracpos;
	  Vector q = wrapfracpos( p );
	  structure.atom.at(i).fracpos = q; //wrapped into cell!
	  //structure.atom.at(i).pos = structure.cell.fractocart( q );
	  smallatom.push_back( structure.atom.at(i) );
  }
  
  //now a huge cell containing many small cells
  
  for ( int superx = p; superx < q+1; superx++ ){
	  for ( int supery = p; supery < q+1 ; supery++ ){
		  for ( int superz = p; superz < q+1 ; superz++ ){
			  for ( int i = 0 ; i < nat ; i++ ){
				  Atom anat = structure.atom.at(i); //copy...
				  //now update the fracpos
				  anat.fracpos.x += superx;
				  anat.fracpos.y += supery;
				  anat.fracpos.z += superz;
				  bigatom.push_back( anat ); //put the atom on the biglist
			  }
		  }
	  }
  }

  //make the Cartesians
  structure.cell.makecell();
  structure.cell.makestar();
  
  for ( int i =0; i < bigatom.size() ; i++) {
    bigatom.at(i).pos = structure.cell.fractocart( bigatom.at(i).fracpos );
    bigatom.at(i).initial_pos = bigatom.at(i).pos;
    cerr << "Atom " << i << " FRAC " << bigatom.at(i).fracpos;
    cerr << " CART " << bigatom.at(i).pos << endl;  
  }
  //all the atoms now have Cartesian pos based on their fractionals
  
  //redefine the basis vectors //using the xi inputs
  Vector oldi = structure.cell.celli;
  Vector oldj = structure.cell.cellj;
  Vector oldk = structure.cell.cellk;
  
  Vector newi = x1*oldi + x2*oldj + x3*oldk;  
  Vector newj = x4*oldi + x5*oldj + x6*oldk;  
  Vector newk = x7*oldi + x8*oldj + x9*oldk;
 
  structure.cell.celli = newi;
  structure.cell.cellj = newj;
  structure.cell.cellk = newk;

  //update the reciprocal vectors to match
  structure.cell.makestar();

  //and update the parameters
  Vector vi = structure.cell.celli;
  Vector vj = structure.cell.cellj;
  Vector vk = structure.cell.cellk;
  //new a
  double newa = sqrt( vi.sq() );
  //new b
  double newb = sqrt( vj.sq() );  
  //new c
  double newc = sqrt( vk.sq() );
  
  //new angles?
  double newcosgamma = vi.dot( vj );
  newcosgamma /= newa;
  newcosgamma /= newb;
  double newgamma = acos( newcosgamma );
  newgamma *= radtodeg; //now in degrees

  double newcosalpha = vj.dot( vk );
  newcosalpha /= newb;
  newcosalpha /= newc;
  double newalpha = acos( newcosalpha );
  newalpha *= radtodeg; //now in degrees

  double newcosbeta = vk.dot( vi );
  newcosbeta /= newa;
  newcosbeta /= newc;
  double newbeta = acos( newcosbeta );
  newbeta *= radtodeg; //now in degrees

  structure.cell.param.at(0) = newa;
  structure.cell.param.at(1) = newb;
  structure.cell.param.at(2) = newc;
  structure.cell.param.at(3) = newalpha;
  structure.cell.param.at(4) = newbeta;
  structure.cell.param.at(5) = newgamma;

  cerr << "New cell parameters: " << endl;
  cerr << structure.cell.printcell() << endl;

  //now, run over bigatom and get new fractional positions from Cartesians!
  for ( int i =0; i < bigatom.size() ; i++) {
    bigatom.at(i).fracpos = structure.cell.carttofrac( bigatom.at(i).pos );
    cerr << "Atom " << i << " FRAC " << bigatom.at(i).fracpos;
    cerr << " CART " << bigatom.at(i).pos << endl;  
  }  
  
  //now, carry out the filter
  vector< Atom > newatom;
  for ( int i = 0 ; i < bigatom.size() ; i++ ){
	  if ( bigatom.at(i).fracpos.x < 0.0 ) continue;
	  if ( bigatom.at(i).fracpos.y < 0.0 ) continue;
	  if ( bigatom.at(i).fracpos.z < 0.0 ) continue;
	  if ( bigatom.at(i).fracpos.x > 1.0 ) continue;
	  if ( bigatom.at(i).fracpos.y > 1.0 ) continue;
	  if ( bigatom.at(i).fracpos.z > 1.0 ) continue;
	  //Vector zz = bigatom.at(i).fracpos;
	  //bigatom.at(i).pos = structure.cell.fractocart( zz );
	  newatom.push_back( bigatom.at(i) );
  }
  
  cerr << "Retaining " << newatom.size() << " atoms in transformed structure." << endl;
  cerr << endl << "CHECK FOR PERIODIC DUPLICATES BEFORE USE!" << endl;
  structure.atom = newatom;

  structure.writecif( outfilename, infilename );

  return(0);
} //ENDOFMAIN
