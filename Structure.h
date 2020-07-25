#include "basic_utils.h"
#include "Vectors.h"
#include "Rotors.h"
#include "Geometry.h"

#ifndef __STRUCTURE_H_INCLUDED_
#define __STRUCTURE_H_INCLUDED_

struct Atom{
  public:

  string species;
  double radius; // keep here for simplicity!
  int role;
  Vector initial_pos;
  Vector pos;
  Vector fracpos;
  vector< int > bondto; // bonded to these atoms
  vector< int > inclust; // in these clusters
  vector< int > isvertex; // vertex id in each cluster

  int type;
  int gridx, gridy, gridz; // coarse grid location

  bool mobile;
  int hybrid; // tag for hybridisation state spN
  //aim to distinguish sp3 from sp2 carbons and nitrogens, for rigidity
  //tag will be 0 in general, 1 for sp2 carbon, 2 for amide nitrogen
  
  vector< int > n14; //list of "1-4" type neighbours for steric scaling
  
  Atom(){
	  species=""; radius = 0.1; role=0; type=0; gridx=0; gridy=0; gridz=0; mobile=true; hybrid=0;
  }
};

const Atom nullatom; // placeholder
bool bonded( Atom &a1, Atom &a2); // checks bonding status using cluster membership

struct Element{
  public:

  string species;
  int role;
  double radius;
  Element(string spin = "", int rolein=0, double radiusin=0.1){
	  species=spin; role=rolein; radius=radiusin;
  }

};

const Element nullelement; // placeholder

class Cell{
	public:
	
	vector< double > param; // the cell parameters a b c alpha beta gamma
	Vector celli, cellj, cellk, istar, jstar, kstar;
	
	double mygridmin; // minimum spacing in coarse grid //set to 3.0 usually
	int Nx, Ny, Nz; // number of cells in coarse grid
	vector< vector< int > > grid; // coarse grid atom list
	
	Cell(){
		mygridmin=3.0; Nx=1; Ny=1; Nz=1;
	}
	
	string printcell(); // output params in string form
	double volume(); // output cell volume using triple product
	void makecell(); // make cell vectors from params
	void makestar(); // make reciprocal vectors
	void makegrid(); // create grid object

	Vector fractocart( Vector& fracin );
	Vector carttofrac( Vector& cartin );
	Vector wrapcartpos( Vector& cartin );
	Vector wrapcartdel( Vector& cartin );

    int posNx( Vector &pos );
    int posNy( Vector &pos );
    int posNz( Vector &pos );
    int fracNx( double fx );
    int fracNy( double fy );
    int fracNz( double fz );
	int Ngrid( int ax, int ay, int az);
	void remove( int index, int gx, int gy, int gz ); //remove index from cell
	void replace( int index, int gx, int gy, int gz ); // put index in cell
		
};

class Structure{
	public:
	
	Cell cell; // all the cell and grid stuff is in here
	vector< Element > element; // my elements live here
	vector< Atom > atom; // my atoms live in here
	vector< Poly > poly;
	vector< Cluster > cluster; // real position clusters
	vector< Cluster > ghost; //ideal position clusters
	
	bool do14; //manage 1-4 neighbours differently in sterics
	double r14scale; //scaling on contact distance for 1-4 neighbours
	
	Structure(){
                do14 = false; r14scale = 1.0;            
    }
	
	bool fillgrid(); // populate the coarse grid with all my atoms
	void formcluster( Cluster &clust ); // make the given cluster's cpos and bonds

	void writextl( string filename, string title, bool do_label = false ); // write out
	void writepol( string filename ); // output polyhedra/clusters
	void writebondlengths( string filename ); // output bond length information
	void writeangles( string filename ); // three-body angles out
	bool readxtl( string filename ); // read in
	bool read_coords( string filename); // read in crystalmaker coord file
	
	vector < int > nearlist( int me, double within ); // atoms near me
	vector< vector < int > > fullnearlist( double within ); // whole-structure near finder	
	Vector mymis( Atom &at , Cluster &clus, int vert );
	void hybridcheck( vector< string > hybridel );
	void mcrigidcheck(); // MC technique
	
	bool elecheck( Atom &atom ); // check an atom against the element list
	bool initialise( double thisgridmin ); // // do initial setup after reading xtl
	void prepclustersandghosts(); // run this to establish geometries of cluster, poly, ghost
	void prepimperfectclustersandghosts(); // run this to establish geometries of cluster, poly, ghost
	void labelvertices(); // run this to label atoms with inclust, isvertex
	void newcell( vector < double > newcellparam ); // update cell, not ghosts
	void newstructure( Structure &news ); // update cell and positions, not ghosts
	
	void fitghosts(); // run fitrotor over the structure's clusters
	double maxrad(); // just get the biggest radius of the elements
	
	vector< Contact > contact( int at , double within );
	//does the Contact list for this atom in this structure using nearlist
	vector< vector< Contact > > fullcontact ( double within );
	//does the Contact list for all atoms, using fullnearlist

	bool bondsintoatoms( vector < Bondline > &bl );
	void polybondfind( vector< Polyspec > polyspec, double pad )	;
	void bondbondfind( vector< Bondspec > bondspec );
	void formbondlines( vector< Bondline > &bl );

	void parsepoly(vector< Polyspec > &polyspec );

	vector< Cluster > candidatecluster ( vector< Bondspec > &bondspec );
	vector< Cluster > Garibaldi( vector< Cluster > &candidate );
	
	void writemismatches( string filename );
	void jiggle( double stepsize ); // random wobble of the structure
	
	void update( vector< Vector > &move ); // update structure using move vVector
	void shift( Vector move ); // shift all atoms by "move" rigidly
	void stripdupes( double prec = 0.0005 ); // removes duplicate atoms;
	void super( int nx, int ny, int nz ); // make supercell
	void othersuper( int nx, int ny, int nz ); // make supercell keeping atom order

    void writecif(string filename, string title); // write out in cif format
    bool readcif( string filename ); //read coords and cell from cif format
    
	void writeclusters( string filename ); // write clusters into an ouput file
	Vector randompos(); // returns a random Cartesian position in the cell
	
	void gridremove( int aid ); // take atom index aid out of cell grid
	void gridreplace( int aid ); // put index aid into cell grid
	int eletype( string &sp ); // give the numeric type of the string species
	
	void find14(); //populate the n14 arrays in the atoms
	bool isn14( int id1, int id2 );
	
};



#endif







