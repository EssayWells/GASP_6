
#include "basic_utils.h"
#include "Vectors.h"
#include "Rotors.h"
#include "Geometry.h"
#include "Structure.h"

#ifndef __GASPCONTROL_H_DEFINED_
#define __GASPCONTROL_H_DEFINED_



extern string inpfilename;
extern bool given_title ;
extern bool hastitle ;
extern string title;
extern bool given_elements ;
extern bool given_poly ;
extern vector< Polyspec > polyspec; // spec when given
extern bool given_bond ;
extern vector< Bondspec > bondspec; // spec when given
extern bool given_option ;
extern bool given_contact ;
extern bool given_grid ;
extern double gridmin; // grid boxes at least this size.
extern double regrid_after; //max move before implementing a regrid
extern double polypad; // multiplier on bondlength when finding poly bonding
extern bool do_jiggle ;
extern double jigsize; //jiggle in Angstroms
extern bool given_relax_command ;
extern bool chebyshev ; // not by default - explodes on contacts! // replace these with int docheb status
extern double cheb; // over-relaxation factor
extern bool use_auto_chebyshev ; //dynamic chebyshev updater
extern double dampclash; // stabilise clash calculations in move
extern bool criterion_mis ;
extern double smallmovecriterion; // stop when move is this small
extern double smallmiscriterion; // and when mis is this small?
extern bool given_input ;
extern bool given_structure ;
extern string input_structure_name;
extern bool given_new_structure ;
extern string new_structure_name;
extern bool given_bond_input ;
extern string input_bonding_name;
extern vector< Bondline > bondline; //bonds from input file if given
extern bool given_new_cell ;
extern bool do_label ; // output labels in the xtl output!
extern bool given_output ;
extern bool do_output_structure ;
extern string output_structure_name;
extern bool do_output_bonding ;
extern string output_bonding_name;
extern bool do_output_poly ;
extern string output_poly_name;
extern bool do_output_bonds ;
extern string output_bonds_name;
extern bool do_output_angles ;
extern string output_angles_name;
extern bool do_output_mismatches ;
extern string output_mismatches_name;
extern bool nodeloc ; // do not include delocalisation rigidity
extern bool use_ion;// special ionic contact rules

extern Structure inputstructure; // use for initial reading
extern Structure structure; // used when actually running
extern Structure newstructure; // used if a restart position is given
extern vector< double > newcellparam; // used if new cell is given

extern bool do_window; //use automatic window search
extern bool found_window; // check before output?
extern string windowtype; // window shape label
extern bool do_outputwindow; //window output to file?
extern string windowbasename; // name for window text and structure files
extern double windowstep; // base step size for search. Set to 0.01 in control.cpp

extern bool do_r14; //modify steric radii for 1-4 cases
extern double r14times; // factor to modify r14 by

extern bool imperfect; //if true, do not make polyhedral geometry ideal

extern bool useXTL; //use XTL instead of (new default) CIF
extern bool useSimple; //simple move calculation instead of the complicated version

extern bool useGradual; //do a gradual move from input to new cell
extern int nGradual; //number of steps
extern string gradualName; //name for gradual output


struct Ioncase{
	public:
	
	string sp1;
	string sp2;
	int type1; // numeric species
	int type2;
	double rclose;
	Ioncase(string sp1in="", string sp2in="", int t1in = -1, int t2in = -2, double rin=1.0){
		sp1 = sp1in, sp2 = sp2in, type1 = t1in, type2 = t2in, rclose = rin;	
	}
	
};
extern vector< Ioncase > ioncase;


//void declareglobal(); // declare all my global variables because C++ hates me
Vector gaspmove( vector< Vector > &mis , vector< Contact > &con, double damp );
Vector simplegaspmove( vector< Vector > &mis , vector< Contact > &con, double damp );
bool gasprelax( Structure &structure, bool &window, double regrid_after=0.1, double damp=0.5, int docheb=0, double cheb=1.0 );
bool validation(); // validates input logic before attempting to run
bool readcommands( string commandfilename); // read the GASP config file

string tokenfromvar( int var );
int dirfromvar( int var );

string randomtag(); //return a tagline from the random list :)

#endif
