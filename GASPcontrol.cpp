#include "GASPcontrol.h"

//GLOBAL VARIABLES instantiated here, listed as extern in header. I hate C++ already.
string inpfilename;
bool given_title = false;
bool hastitle = false;
string title;
bool given_elements = false;
bool given_poly = false;
vector< Polyspec > polyspec; // spec when given
bool given_bond = false;
vector< Bondspec > bondspec; // spec when given
bool given_option = false;
bool given_contact = false;
bool given_grid = false;
double gridmin = 3.0; // grid boxes at least this size.
//gridmin = 3.0;
double regrid_after = 0.1; //max move before implementing a regrid
double polypad = 1.3; // multiplier on bondlength when finding poly bonding
bool do_jiggle = false;
double jigsize = 0.1; //jiggle in Angstroms //implement!
bool given_relax_command = false;
//note chebyshev now deprecated and not activated
//remove in next revision
bool chebyshev = false; // not by default - explodes on contacts! // replace these with int docheb status
double cheb = 1.0; // over-relaxation factor
bool use_auto_chebyshev = false; //dynamic chebyshev updater
double dampclash = 0.5; // stabilise clash calculations in move
bool criterion_mis = false;
double smallmovecriterion = 1E-6; // stop when move is this small
double smallmiscriterion = 1E-3; // and when mis is this small?
bool given_input = false;
bool given_structure = false;
string input_structure_name;
bool given_new_structure = false;
string new_structure_name;
bool given_bond_input = false;
string input_bonding_name;
vector< Bondline > bondline; //bonds from input file if given
bool given_new_cell = false;
bool do_label = false; // output labels in the xtl output!
bool given_output = false;
bool do_output_structure = false;
string output_structure_name;
bool do_output_bonding = false;
string output_bonding_name;
bool do_output_poly = false;
string output_poly_name;
bool do_output_bonds = false;
string output_bonds_name;
bool do_output_angles = false;
string output_angles_name;
bool do_output_mismatches = false;
string output_mismatches_name;
bool nodeloc = false; // do not include delocalisation rigidity

Structure inputstructure; // use for initial reading
Structure structure; // used when actually running
Structure newstructure; // used if a restart position is given
vector< double > newcellparam; // used if new cell is given

bool do_window = false; //use automatic window search
bool found_window = false;
string windowtype; // window shape label
bool do_outputwindow = false; //window output to file?
string windowbasename; // name for window text and structure files
double windowstep = 0.01; // base step size for search.

bool do_r14 = false ; //modify steric radii for 1-4 cases
double r14times = 1.0 ; // factor to modify r14 by

bool imperfect = false; // by default, do not make polyhedra imperfect
bool useXTL = false; //use XTL instead of (new default) CIF
bool useSimple = false; // use simple move calculation

vector< Ioncase > ioncase;
bool use_ion = false; // use special ionic rules in contact

bool useGradual = false; //do a gradual move from input to new cell
int nGradual = 1; //number of steps
string gradualName = "gradual" ; //name for gradual output



bool validate_ioncase( Structure &structure, Ioncase &ionc ){
	//check string species versus elements
	//assign numeric type
	//return false on unrecognised element
	
	ionc.type1 = structure.eletype( ionc.sp1 );
	ionc.type2 = structure.eletype( ionc.sp2 );
	
	if ( ionc.type1 < 0 || ionc.type2 < 0 ) return false;
	return true;
}

//
//gaspmove considers the vertex and contact positions
//and iterates to a new "relaxed" position for the central atom
//result is a move to be applied to the atom when we update
//
Vector gaspmove( vector< Vector > &mis , vector< Contact > &con, double damp=0.5 ){
       Vector xold = Vector(0,0,0); // old x
       Vector xnew = Vector(0,0,0); // new x
       
       Vector xorig = Vector(0,0,0); // use to soften contacts
       
       //double smalld2 = 1e-16; // when difference becomes "small"
       double smalld2 = 1e-16; // when difference becomes "small"
       double xn, yn, zn; // track components;
       
       int nc = con.size();
       int nm = mis.size();
       
       if ( nc == 0 && nm == 0 ){
            return xold; //nothing to be done!
       }

       if ( nm > 0 && nc > 0 ){
          for ( int m = 0 ; m < nm ; m++ ){
              xold += mis.at(m);
          }
          xold /= nm; // preset to vertex-matching position before starting clash calculation
       }
       xorig = xold;
       
       bool working = true;
       int counter = 0;
       int maxcount = 100;
       
       while (working){
             counter++;
             //cerr << "Iteration " << counter << ": ";
             double xu, xd; // numerator and denominator
             double yu, yd;
             double zu,zd;

             xd = yd = zd = nm;
             xu = yu = zu = 0.0;
             for ( int m = 0 ; m < nm ; m++ ){
                 xu +=mis.at(m).x;
                 yu +=mis.at(m).y;
                 zu +=mis.at(m).z;
             }

             // do dfactors
             double dfac;
             double dist;
             double ratio;
             Vector diff;
             for ( int c = 0; c < nc ; c++ ){
                 diff = con.at(c).dr - xold; // current vector
                 dist = sqrt( diff.sq() );
                 ratio = con.at(c).l / dist ;
                 //if ( ratio > 1.1 ){
                 //     cerr << "WARNING: clash worse than 1.1." << endl;
                 //     ratio = 1.1; // trap?
                 //}
                 dfac = 1 - ratio ; //weighting on contact
                 xd += dfac * damp;
                 yd += dfac * damp;
                 zd += dfac * damp;
                 //xu += con.at(c).dr.x * dfac;
                 //yu += con.at(c).dr.y * dfac;
                 //zu += con.at(c).dr.z * dfac;
                 xu += con.at(c).dr.x * dfac * damp;
                 yu += con.at(c).dr.y * dfac * damp;
                 zu += con.at(c).dr.z * dfac * damp;
             }
             
             //new approach: trap if contacts exist
             //effectively starting position is a pseudovertex
             if ( nc > 0 ){
                  xu += xorig.x;
                  yu += xorig.y;
                  zu += xorig.z;
                  xd += 1;
                  yd += 1;
                  zd += 1;
             }
             
             //new x
             //cerr << "UP: " << xu << " " << yu << " " << zu << endl;
             //cerr << "DN: " << xd << " " << yd << " " << zd << endl;
             xn = xu / xd;
             yn = yu / yd;
             zn = zu / zd;
             xnew = Vector ( xn, yn, zn );
             //cerr << "Xnew = " << xnew << endl;
             diff = xnew - xold;
             
             if ( diff.sq() < smalld2 ) working = false;
             if ( counter > maxcount ){
                  working = false; // bored now
                  //cerr << "WARNING: move calculation not well converged." << endl;
                  //xnew = Vector(0,0,0); // don't send a badly converged move
             }
             xold = xnew;
       }
       
       return xold;
}

//
//simplegaspmove considers the vertex and contact positions
//and averages a response position for the atom with no cleverness
//result is a move to be applied to the atom when we update
//
Vector simplegaspmove( vector< Vector > &mis , vector< Contact > &con, double damp=0.5 ){
	   //damp is not used!
       Vector xold = Vector(0,0,0); // old x = here
       Vector xnew = Vector(0,0,0); // new x
       
       int nc = con.size();
       int nm = mis.size();
       int ntotal = nc + nm;
       
       if ( nc == 0 && nm == 0 ){
            return xold; //nothing to be done!
       }
       
       //if I have a contact but no bonds (unlikely?)
       //then pass a fictitious brake
       //by counting a zero vector into xnew
       //as a fake constraint
       if ( nm == 0 && nc > 0 ) ntotal++; 

       for ( int m = 0 ; m < nm ; m++ ){
           xnew += mis.at(m);
       }
       //we have counted each mismatch vector as a contribution
       //now loop over contacts obtaining their vector
       for ( int c = 0; c < nc ; c++ ){
                 Vector diff = con.at(c).dr - xold; // vector to contacting atom
                 //how long is it?
                 double dist = sqrt( diff.sq() );
                 //so how long is the overlap?
                 double overlap = con.at(c).l - dist;
                 //so my vector contribution
                 //is antiparallel to diff
                 //and sized as overlap
                 diff.norm(); //unit sized now
                 diff *= overlap;
                 xnew -= diff; //vector component done
       }
       
       //an average of everything:
       xnew /= ntotal;
       
       return xnew;
}


//
//gasprelax makes the structure do the GASP relaxation loop
//smallmove and smallmis criteria are provided by GASP variables, or should be!
//
 bool gasprelax( Structure &structure, bool &window,  double regrid_after, double damp, int docheb, double cheb  ){
    window = false; // criterion check
    bool usetrend = false; //use for testing "trend" system without breaking everything
    usetrend = true; //use for testing "trend" system without breaking everything
	vector< Vector > lastmove, move;
	vector< Vector > trend;
	move.resize( structure.atom.size() );
	if ( docheb > 1 ) lastmove.resize( structure.atom.size() ); // use for tracking if cheb auto
	trend = move;
	
	cerr << "Relaxing: " << endl << endl;
	
	bool fitting = true;
	bool becausemove = false;
	bool becausemis = false; // move or mistmatch convergence
	
	double movetrack = 0.0;
	double m2 = 0.;
	double bigm2 = 0.;
	double bigm = 0.;
	
	double worstm = 0.; // worst distortion
	double worstc = 0.; // worst steric
	
	int counter=0;
	//fairly tight scaling for starters!
	int startscaling = 100; //start to damp oscillations if we got this far
    int endscaling = 100000;
	int whichbigmis=0, whichbigc=0, whichbigmove=0, whichbigmisclus=0, whichbigmisvert=0;	
	
	while( fitting ){
		counter++; 
		cerr << "Fitting cycle " << fixed << setw(8) << counter << "; ";
		bigm2 = 0.;
		
		structure.fitghosts(); // ghosts rotate to fit cluster
		
		worstm=0.;
		worstc=0.;
		double lookie = 2.0*structure.maxrad();
		
		//1-4 issues will be dealt with implicitly inside contact
		vector< vector< Contact > > thecon = structure.fullcontact( lookie );

		
		for ( int i = 0; i < structure.atom.size(); i++){
			vector< Vector > miss;
			miss.clear();
			Vector amis;
			int ncl = structure.atom.at(i).inclust.size();
			
			//mismatches
			for ( int wc = 0; wc < ncl; wc++ ){
				int ac = structure.atom.at(i).inclust.at(wc);
				int av = structure.atom.at(i).isvertex.at(wc);
				amis = structure.mymis( structure.atom.at(i), structure.ghost.at(ac), av );
				miss.push_back( amis );
				
				double dm = amis.sq();
				if ( dm > worstm ){
					worstm = dm;
					whichbigmis = i;
					whichbigmisclus = ac;
					whichbigmisvert = av;
				}
			}
			
			//contacts
			vector< Contact > cont = thecon.at(i);
			int ncon = cont.size();
			for ( int wcon = 0; wcon < ncon; wcon++ ){
				double c = cont.at(wcon).l - sqrt( cont.at(wcon).dr.sq() );
				if ( c > worstc ) worstc = c;
				if ( c > worstc ) whichbigc = i;
			}
			
			//if using ion option, mess with the collisions here
			if ( use_ion ){
				for ( int wcon = 0; wcon < ncon; wcon++ ){
					int a1 = cont.at(wcon).id1;
					int a2 = cont.at(wcon).id2;
					int t1 = structure.atom.at(a1).type;
					int t2 = structure.atom.at(a2).type;
					for ( int anion = 0; anion < ioncase.size(); anion++ ){
						int ot1 = ioncase.at(anion).type1;
						int ot2 = ioncase.at(anion).type2;
						if ( ( ot1 == t1 && ot2 == t2 ) || ( ot1 == t2 && ot2 == t1 ) ){
							//reset this contact distance
							cont.at(wcon).l = ioncase.at(anion).rclose;
							break;
						}
					}
				}
			}
			
			Vector amove;
			if ( useSimple ){
				amove = simplegaspmove( miss, cont, damp );
			}
			else{
				amove = gaspmove ( miss, cont, damp ); // explicit damping parameter
			}
			
			//okay now implement emergency damping function
			//scale down the moves if the counter is getting too high
	
			if ( counter > startscaling ){
				double factor = ( counter - startscaling ) / endscaling;
				double scaler = 1. - factor;
				amove *= scaler;
				//so amove will start to decline gently from 1K to zero at 100K
				//so smallmove will break the loop at some point
				//should guarantee stability at least
			}
			
			move.at(i) = amove; // into array
			m2 = amove.sq();
			if ( m2 > bigm2 ){
				bigm2 = m2; whichbigmove = i;
			}
			
			if ( usetrend){
				trend.at(i) *= 0.95;
				//trend.at(i) += 0.1*amove;
				trend.at(i) += amove;
				//trend is a "running" count of how the moves are adding up over time
			}
			
		}
		bigm = sqrt( bigm2 );
		double bigmis = sqrt( worstm );
		//if we need to trap com move, do it here
        cerr << "Mismatch; " << fixed << setw(12) << setprecision(7) << bigmis;
        cerr << " ; Clash " << fixed << setw(12) << setprecision(7) << worstc << " ; ";
        cerr << "Move: " << fixed << setw(14) << setprecision(10) << bigm << endl;
        /*
        if ( usetrend ){
			cerr << "MOVE " << move.at( whichbigmove );
			cerr << " TREND " << trend.at( whichbigmove );
			cerr << " FOR " << whichbigmove+1 << endl;
		}
		*/
		
		//criteria
        bool stopmove = false;
        bool stopmis = false;
             
        if ( !criterion_mis ){
			stopmis = false; // no mis criterion, ignore
        }
		else{
			if ( bigmis < smallmiscriterion ) stopmis = true;
		}

		if ( bigm < smallmovecriterion ) stopmove = true;

		becausemove = stopmove;
		becausemis = stopmis;

		if ( stopmis || stopmove ) fitting = false; // done fitting
		if ( counter >= endscaling ) fitting = false; //I'm bored now 		
		
		//cheb
		if ( cheb > 0  && fitting && counter > 10 ){
			if ( cheb == 1 ){
				bigm *= cheb; //scaling
				for ( int i = 0; i < structure.atom.size() ; i++){
					move.at(i) *= cheb; // scaling
				}
			}
			else if ( cheb == 2 ){ // cheb auto logic
				double totdot = 0.;
				double thisdot = 0.;
				int npos = 0;
				int nneg = 0;
				for ( int i = 0 ; i < structure.atom.size(); i++ ){
					thisdot = move.at(i).dot( lastmove.at(i) );
					totdot += thisdot;
					if ( thisdot < 0. ){ nneg ++ ; }
					else { npos++ ; }
				}  				

                if ( thisdot < 0. ){
                     //system is in oscillation: damp it
                     cheb = 0.5;
                }
                else if ( nneg > 0 ){
                     //some atoms are oscillating
                     //run default
                     cheb = 1.0;
                }
                else if ( bigm > 0.01 ){
                     //the step size is quite large
                     cheb = 1.0;
                }
                else {
                     //accelerate! Ramming speed!
                     cheb = 1.9;
                }
				bigm *= cheb; //scaling
				for ( int i = 0; i < structure.atom.size() ; i++){
					move.at(i) *= cheb; // scaling
				}
				lastmove = move; // cycle for next time				
			}
		}

		if ( usetrend && counter > 1000 ){
			structure.update( trend );
		}
		else{
			structure.update( move );
		}
		movetrack += bigm; // cumulative largest step

		if ( movetrack > regrid_after){
			movetrack = 0.0;
			cerr << endl;
			cerr << "Regridding." << endl;
			structure.cell.makegrid();
			bool goodgrid = false;
			goodgrid = structure.fillgrid();
			if (!goodgrid){
				cerr << "FATAL PROBLEM: grid failure." << endl;
				window = false;
				return false; // we failed to relax
			}
		}		
		
	} // finished fitting
	cerr << "Fitting complete:" << endl;
	if ( becausemove ) cerr << "Move less than criterion " << smallmovecriterion << endl;
	if ( becausemis ) cerr << "Mismatch less than criterion " << smallmiscriterion << endl;            
	if ( !becausemove && ! becausemis ) cerr << "Because I said so." << endl;	
	
	if (becausemis) window = true; // in a window.
	
	return true; // we finished the relaxation
}

//
//validation routine checks for a valid input before running
//
bool validation(){
     cerr << "Entering validation process." << endl;
     if ( ! given_input ) {
        cerr << "FATAL: No input block given." << endl;
        return false;
     }     
     if ( ! given_structure ) {
        cerr << "FATAL: No input structure given." << endl;
        return false;
     }     
     if ( ! given_elements || ( inputstructure.element.size() == 0 ) ) {
        cerr << "FATAL: No elements given." << endl;
        cerr << "Leeloo Dallas multipass." << endl;
        return false;
     }
     if ( ! given_output ) {
        cerr << "WARNING: no output requested. Really?" << endl;
     }    
     if ( ( ! given_bond) && ( ! given_poly ) ) {
        cerr << "FATAL: No cluster-building information." << endl;
        cerr << "Please provide at least one of: BOND block, POLY block." << endl;
        return false;
     }

     if ( ( given_new_structure || given_new_cell ) && ! given_relax_command ) {
        cerr << "FATAL: must request relax to adapt to new structure or cell." << endl;
        return false;     
     }
     if ( given_bond && ( bondspec.size() == 0 ) ){
          cerr << "WARNING: BOND block given but I have no bond specifications. Really?" << endl;
     }
     if ( given_poly && ( polyspec.size() == 0 ) ){
          cerr << "WARNING: POLY block given but I have no poly specifications. Really?" << endl;
     }
     if ( bondspec.size() == 0 && polyspec.size() == 0 ) {
          cerr << "FATAL: no bond or poly specifications given, cannot build clusters!" << endl;
          return false;
     }   
     if ( use_ion ){
		 cerr << "Using ION option: validating ion cases." << endl;
		 cerr << "Got " << ioncase.size() << " ion cases." << endl;
		 for ( int i = 0 ; i < ioncase.size(); i++ ){
			 bool anygood = validate_ioncase( inputstructure, ioncase.at(i) );
			 if (!anygood ){
				 cerr << "FATAL: unrecognised element type in ioncase " << i << endl;
				 cerr << "   " << ioncase.at(i).sp1 << "," << ioncase.at(i).sp2 << endl;
				 return false;
			 }
		 }
	 }
	 if (do_window && !given_relax_command){
		 cerr << "WARNING: Window search option implies relax option: activating." << endl;
		 given_relax_command = true;
	 }
	 if ( useGradual && !given_relax_command){
		 cerr << "WARNING: GRADUAL option implies relax option: activating." << endl;
		 given_relax_command = true;
	 }
	 if ( useGradual && do_window ){
		 cerr << "ERROR: GRADUAL option is incompatible with SEARCH option." << endl;
		 return false;
	 }
	 if ( do_window && !do_outputwindow ){
		 cerr << "WARNING: asked for window search, but no window file ouput requested. Really?" << endl;
	 }
	 if ( !do_window && do_outputwindow ){
		 cerr << "ERROR: can't ask for window file output without window search option." << endl;
		 return false;
	 }
	 if ( do_window && !criterion_mis ){
		 cerr << "ERROR: window search requires SMALLMISmatch criterion." << endl;
		 return false;
	 }

 cerr << "Input logic valid, proceeding." << endl;
 return true; // passed logic tests    
}

//
//read the config file for GASP
//
bool readcommands( string commandfilename){
     string linein;
     string foo; //use as buffer
     vector< string > tokensin, Utokensin; //U is not lowercased :)
     vector< string > unhashed, Uunhashed; //Yes Uunhashed
     ifstream inputfile( commandfilename.c_str() );
     
     if ( inputfile.fail() ) {
        cerr << "Problem! Cannot read file " << commandfilename << endl;
        return false;     
     }

     bool isinblock = false; // are we in a defined block?
     //keywords: title element input bond poly option output
     bool intitle = false;
     bool inelement = false;
     bool ininput = false;
     bool inbond = false;
     bool inpoly = false;
     bool inoption = false;
     bool inoutput = false;
     
     bool iskey = false; // spot keywords
     bool isslash = false; // spot block-ending slash markers
     double bondwithin = 2.0; //check for labelled bond blocks!
     double defaultbondwithin = 2.0;

     while ( getline(inputfile, foo) ){
           iskey = false;
           isslash = false;
           
           string woof = foo;
           linein = to_lower( woof ) ; // all lower case for simplicity
           //but we also make Utokensin from foo for filename nonlowercasing!
           //cerr << endl;
           //cerr << foo << endl;
           //cerr << linein << endl;
           
           
           tokensin = chop(linein, " \t"); // array of tokens
           Utokensin = chop(foo, " \t"); // array of tokens from foo
           unhashed.clear();
           Uunhashed.clear();
           if ( tokensin.size() < 1 ){
              continue; // blank line
           }
           // drop everything after the first hash token
           for ( int i = 0 ; i < tokensin.size() ; i++ ) {
               string hasher = tokensin.at(i).substr(0,1); // first character
               if ( hasher == "#" ) break; // leave for loop now
               unhashed.push_back( tokensin.at(i) ); // unhashed skips all hash
               Uunhashed.push_back( Utokensin.at(i) ); // non-lowercase version for filenames
           }
           if ( unhashed.size() < 1 ){
              continue; // nothing left but hashes and whitespace
           }
           //if we're here, we have something to do

           if ( isinblock ){
              //check for block-ending tag
              string wordone = unhashed.at(0);
              string slasher = wordone.substr(0,1); // first character of first word
              if ( slasher == "/" ) {
                 isslash = true; // now what?
                 //cerr << "Slasher movie at line: " << foo << endl;         
                 //Check if this is the element ending the current block, or a mistake
                 if ( wordone == "/title"){
                      if ( intitle ) {
                           isinblock = false;
                           intitle = false;
                           continue; // go for next line        
                      }
                      else {
                           //wrong block!
                           cerr << "FATAL: /title keyword, not in title block." << endl;
                           return false;
                      }
                 }
                 if ( wordone == "/element"){
                      if ( inelement ) {
                           isinblock = false;
                           inelement = false;
                           continue; // go for next line        
                      }
                      else {
                           //wrong block!
                           cerr << "FATAL: /element keyword, not in element block." << endl;
                           return false;
                      }
                 }
                 if ( wordone == "/input"){
                      if ( ininput ) {
                           isinblock = false;
                           ininput = false;
                           continue; // go for next line        
                      }
                      else {
                           //wrong block!
                           cerr << "FATAL: /input keyword, not in input block." << endl;
                           return false;
                      }
                 }
                 if ( wordone == "/bond"){
                      if ( inbond ) {
                           isinblock = false;
                           inbond = false;
                           continue; // go for next line        
                      }
                      else {
                           //wrong block!
                           cerr << "FATAL: /bond keyword, not in bond block." << endl;
                           return false;
                      }
                 }
                 if ( wordone == "/poly"){
                      if ( inpoly ) {
                           isinblock = false;
                           inpoly = false;
                           continue; // go for next line        
                      }
                      else {
                           //wrong block!
                           cerr << "FATAL: /poly keyword, not in poly block." << endl;
                           return false;
                      }
                 }                 
                 if ( wordone == "/option"){
                      if ( inoption ) {
                           isinblock = false;
                           inoption = false;
                           continue; // go for next line        
                      }
                      else {
                           //wrong block!
                           cerr << "FATAL: /option keyword, not in option block." << endl;
                           return false;
                      }
                 }                    
                 if ( wordone == "/output"){
                      if ( inoutput ) {
                           isinblock = false;
                           inoutput = false;
                           continue; // go for next line        
                      }
                      else {
                           //wrong block!
                           cerr << "FATAL: /output keyword, not in output block." << endl;
                           return false;
                      }
                 }                 
              } // that's the end of the slash check
              //check for data according to block type
              
              if ( intitle ) {
                   //title processing
                   title = linein;
                   hastitle = true;
                   cerr << "Title: " << linein << endl;
              }
              else if ( inelement ) {
                   //element processing
                   if ( unhashed.size() > 1 ) {
                        char dig = unhashed.at(1).at(0);
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in element block." << endl;
                             continue; // not suitable to be an element spec     
                        }
                        Element newelement;
                        newelement.species = unhashed.at(0);
                        newelement.role = 0;
                        newelement.radius = atof( unhashed.at(1).c_str() );
                        inputstructure.element.push_back( newelement );
                        //if ( newelement.radius > maxrad ) maxrad = newelement.radius; // track biggest sphere
                        cerr << "Found element " << newelement.species << " , " << newelement.radius << endl;
                   }
              }   
              else if ( ininput ) {
                   //input processing
                   if ( unhashed.size() < 2 ){
                        cerr << "Skipping line " << linein << " in input block." << endl;
                        continue; // all entries must be two or more tokens!
                   }
                   if ( unhashed.at(0) == "structure" ) {
                        given_structure = true;
                        input_structure_name = Uunhashed.at(1);
                        cerr << "Will read input structure from " << input_structure_name << endl;
                        continue; // read next line
                   }
                   if ( unhashed.at(0) == "bonding" ) {
                        given_bond_input = true;
                        input_bonding_name = Uunhashed.at(1);
                        cerr << "Will read input bonding from " << input_bonding_name << endl;
                        continue;
                   }    
                   if ( unhashed.at(0) == "new" ){
                        if ( unhashed.at(1) == "structure" ) {
                             if ( unhashed.size() > 2 ) {
                                  given_new_structure = true;
                                  new_structure_name = Uunhashed.at(2);
                                  cerr << "New structure given in file " << new_structure_name << endl;
                                  continue;
                             }
                             else {
                                  cerr << "Skipping line " << linein << " in input block." << endl;
                                  continue;
                             }
                        }
                        else if ( unhashed.at(1) == "cell" ) {
                             if ( unhashed.size() < 8 ) {
                                  cerr << "Skipping line " << linein << " in input block." << endl;
                                  continue;
                             }
                             else{
                                  cerr << "New cell params from line " << linein << endl;
                                  given_new_cell = true;
                                  newcellparam.clear();
                                  for ( int i = 0 ; i < 6 ; i++ ) {
                                      double param = atof( unhashed.at(i+2).c_str() );
                                      newcellparam.push_back( param );
                                  }
                                  continue;
                             }
                        }
                   }               
                   
              }           
              else if ( inbond ) {
                   //bond processing
                   if ( unhashed.size() > 1 ) {
                        cerr << "Bonding: " << linein << endl;
                        Bondspec newbondspec;
                        newbondspec.within = bondwithin;
                        newbondspec.speciesA = unhashed.at(0);
                        newbondspec.speciesB.clear();
                        for ( int i = 1 ; i < unhashed.size() ; i++ ){
                            newbondspec.speciesB.push_back( unhashed.at(i) ) ;
                        }
                        bondspec.push_back( newbondspec );
                   }
              }               
              else if ( inpoly ) {
                   //poly processing
                   if ( unhashed.size() > 3 ) {
                        char dig = unhashed.at(3).at(0); // fourth item should be a bondlength!
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in poly block." << endl;
                             continue; // not suitable to be a poly spec     
                        }
                        string key = unhashed.at(0).substr(0,3);
                        if ( key == "tet" ){
                             //tetrahedron
                             Polyspec newpolyspec;
                             newpolyspec.shape = 1; // type 1 is tetrahedron
                             newpolyspec.c_species = unhashed.at(1);
                             newpolyspec.v_species = unhashed.at(2);
                             newpolyspec.bondlength = atof( unhashed.at(3).c_str() );
                             polyspec.push_back ( newpolyspec);
                             cerr << "Found TETrahedron on line: " << linein << endl;
                             continue;
                        }
                        if ( key == "squ" ){
                             //square
                             Polyspec newpolyspec;
                             newpolyspec.shape = 4; // type 4 is square planar
                             newpolyspec.c_species = unhashed.at(1);
                             newpolyspec.v_species = unhashed.at(2);
                             newpolyspec.bondlength = atof( unhashed.at(3).c_str() );
                             polyspec.push_back ( newpolyspec);
                             cerr << "Found SQUare on line: " << linein << endl;
                             continue;
                        }
                        if ( key == "oct" ){
                             //octahedron
                             Polyspec newpolyspec;
                             newpolyspec.shape = 2; // type 2 is octahedron
                             newpolyspec.c_species = unhashed.at(1);
                             newpolyspec.v_species = unhashed.at(2);
                             newpolyspec.bondlength = atof( unhashed.at(3).c_str() );
                             polyspec.push_back ( newpolyspec);
                             cerr << "Found OCTrahedron on line: " << linein << endl;
                             continue;
                        }
                        if ( key == "tri" ){
                             //triangle
                             Polyspec newpolyspec;
                             newpolyspec.shape = 3; // type 3 is triangle
                             newpolyspec.c_species = unhashed.at(1);
                             newpolyspec.v_species = unhashed.at(2);
                             newpolyspec.bondlength = atof( unhashed.at(3).c_str() );
                             polyspec.push_back ( newpolyspec);
                             cerr << "Found TRIangle on line: " << linein << endl;
                             continue;
                        }                        
                        if ( key == "bar" ){
                             //bar
                             Polyspec newpolyspec;
                             newpolyspec.shape = 5; // type 5 is bar (diatomic)
                             newpolyspec.c_species = unhashed.at(1);
                             newpolyspec.v_species = unhashed.at(2);
                             newpolyspec.bondlength = atof( unhashed.at(3).c_str() );
                             polyspec.push_back ( newpolyspec);
                             cerr << "Found BAR on line: " << linein << endl;
                             continue;
                        }                        
                   }
              }               
              else if ( inoption ) {
                   //option processing
                   //contact; relax; others?
                   if ( unhashed.at(0) == "relax" ) {
                        cerr << "RELAX option selected." << endl;
                        given_relax_command = true;
                        continue;
                   }
                   if ( unhashed.at(0) == "xtl" ) {
                        cerr << "Use XTL option selected." << endl;
                        useXTL = true;
                        continue;
                   }                   
                   if ( unhashed.at(0) == "simple" ) {
                        cerr << "SIMPLE move calculation option selected." << endl;
                        useSimple = true;
                        continue;
                   }
                   if ( unhashed.at(0) == "gradual" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             cerr << "GRADUAL option needs an integer argument." << endl;
                             continue; // not suitable     
                        }
                        useGradual = true;
                        nGradual = atoi ( unhashed.at(1).c_str() );
                        cerr << "Setting GRADUAL steps to " << nGradual << endl;
                   }
                   if ( unhashed.at(0) == "imperfect" ) {
                        cerr << "IMPERFECT option selected." << endl;
                        imperfect = true;
                        continue;
                   }                   
                   //any more options, put here!
                   if ( unhashed.at(0) == "nodeloc" ) {
                        cerr << "NODELOC option selected; will not unify sp2 centers." << endl;
                        nodeloc = true;
                        continue;
                   }
                   if ( unhashed.at(0) == "grid" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        given_grid = true;
                        gridmin = atof ( unhashed.at(1).c_str() );
                        cerr << "Setting GRID box size to " << gridmin << endl;
                   }
                   if ( unhashed.at(0) == "search" || unhashed.at(0) == "window" ){
					   if ( unhashed.size() ==1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
					   }
					   string wintype = unhashed.at(1).substr(0,3); //should label a shape type
					   if ( wintype == "cub" || wintype == "iso" || wintype == "six" || wintype == "tet" || wintype == "hex" || wintype == "ort" || wintype == "rho" || wintype == "mon" || wintype == "tri" || wintype == "aaa" || wintype == "bbb" || wintype == "ccc" || wintype == "alp" || wintype == "bet" || wintype == "gam" ){
						   windowtype = wintype;
						   do_window = true;
						   cerr << "Using WINDOW SEARCH option with window type " << windowtype << endl;
					   }
					   else{
						   cerr << "ERROR: unrecognised window type " << wintype << endl;
						   cerr << "NOT using WINDOW SEARCH option, sorry." << endl;   
					   }
				   }             
                   if ( unhashed.at(0) == "polypad" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        //given_grid = true;
                        polypad = atof ( unhashed.at(1).c_str() );
                        cerr << "Setting POLYhedron PADding factor to " << polypad << endl;
                   }
                   if ( unhashed.at(0) == "cheb" ) {
					   cerr << "No longer using previous CHEByshev option." << endl;
					   cerr << "cheb option not active." << endl;
					   /*
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        //detect auto chebyshev command
                        if ( unhashed.at(1) == "auto" ){
                             cerr << "Detected AUTO CHEByshev option." << endl;
                             chebyshev = true; 
                             use_auto_chebyshev = true;
                             continue; // done, use auto setting
                        }
                        //if we're still here, we want a number for the argument
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        chebyshev = true;
                        cheb = atof ( unhashed.at(1).c_str() );
                        cerr << "Setting CHEByshev acceleration factor to " << cheb << endl;
                        */ 
                   }
                   if ( unhashed.at(0) == "smallmove" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        //given_grid = true;
                        smallmovecriterion = atof ( unhashed.at(1).c_str() );
                        cerr << "Setting SMALLMOVE criterion to " << smallmovecriterion << endl;
                   }
                   if ( unhashed.at(0) == "smallmis" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        //given_grid = true;
                        smallmiscriterion = atof ( unhashed.at(1).c_str() );
                        criterion_mis = true;
                        cerr << "Activating SMALLMISmatch criterion." << endl;
                        cerr << "Setting SMALLMISmatch criterion to " << smallmiscriterion << endl;
                   }
                   if ( unhashed.at(0) == "dampclash" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        //given_grid = true;
                        dampclash = atof ( unhashed.at(1).c_str() );
                        cerr << "Setting SMALLMISmatch criterion to " << dampclash << endl;
                   }
                   if ( unhashed.at(0) == "jiggle" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        jigsize = atof ( unhashed.at(1).c_str() );
                        do_jiggle = true;
                        cerr << "Setting JIGGLE size to " << jigsize << endl;
                   }
                   if ( unhashed.at(0) == "r14" ) {
                        if ( unhashed.size() == 1 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        char dig = unhashed.at(1).at(0); // second item should be a double
                        if ( ! isdigit( dig ) ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue; // not suitable     
                        }
                        r14times = atof ( unhashed.at(1).c_str() );
                        do_r14 = true;
                        cerr << "Setting R14 scale factor to " << r14times << endl;
                   }                   
                   if ( unhashed.at(0) == "ion" ) {
                        if ( unhashed.size() < 4 ){
                             cerr << "Skipping line " << linein << " in options block." << endl;
                             continue;
                        }
                        use_ion = true;
                        Ioncase ic;
                        string s1 = unhashed.at(1).c_str();
                        string s2 = unhashed.at(2).c_str();
                        double r = atof( unhashed.at(3).c_str() );
                        ic.sp1 = s1;
                        ic.sp2 = s2;
                        ic.rclose = r;
                        ioncase.push_back( ic );
                        cerr << "Setting ION option on " << s1 << "," << s2 << "," << r << endl;
                   }                   
                   if ( unhashed.at(0) == "label" ) {
                        do_label = true;
                        cerr << "LABEL on: will put atom IDs after coords in xtl output." << jigsize << endl;
                   }                   
              }               
              else if ( inoutput ) {
                   //output processing
                   if ( unhashed.size() < 2 ){
                        cerr << "Skipping line " << linein << " in output block." << endl;
                        continue; // all entries must be two tokens!
                   }
                   if ( unhashed.at(0) == "structure" ) {
                        do_output_structure = true;
                        output_structure_name = Uunhashed.at(1);
                        cerr << "Will write output structure to " << output_structure_name << endl;
                        continue; // read next line
                   }
                   if ( unhashed.at(0) == "window" ) {
                        do_outputwindow = true;
                        windowbasename = Uunhashed.at(1);
                        cerr << "Will write window edge structures with base name: " << windowbasename << endl;
                        continue; // read next line
                   }
                   if ( unhashed.at(0) == "gradualname" ) {
                        gradualName = Uunhashed.at(1);
                        cerr << "Will write GRADUAL structures with base name: " << gradualName << endl;
                        continue; // read next line
                   }                    
                   if ( unhashed.at(0) == "bonding" ) {
                        do_output_bonding = true;
                        output_bonding_name = Uunhashed.at(1);
                        cerr << "Will write output bonding topology to " << output_bonding_name << endl;
                        continue;
                   }                   
                   if ( unhashed.at(0) == "pol" ) {
                        do_output_poly = true;
                        output_poly_name = Uunhashed.at(1);
                        cerr << "Will write output polyhedra to " << output_poly_name << endl;
                        continue;
                   }           
                   if ( unhashed.at(0) == "bonds" || unhashed.at(0) == "bondlengths" ) {
                        do_output_bonds = true;
                        output_bonds_name = Uunhashed.at(1);
                        cerr << "Will write bond lengths to " << output_bonds_name << endl;
                        continue;
                   } 
                   if ( unhashed.at(0) == "angles" ) {
                        do_output_angles = true;
                        output_angles_name = Uunhashed.at(1);
                        cerr << "Will write bond lengths to " << output_angles_name << endl;
                        continue;
                   } 
                   if ( unhashed.at(0) == "mismatch" || unhashed.at(0) == "clash" ) {
                        do_output_mismatches = true;
                        output_mismatches_name = Uunhashed.at(1);
                        cerr << "Will write MISMATCH and CLASH info to " << output_mismatches_name << endl;
                        continue;
                   }
              }           
           }
           else {
                // not in a block
                //check for block starting tag, ignore everything else
                string wordone = unhashed.at(0);
                if ( wordone == "title" ) {
                     isinblock = true;
                     intitle = true;
                     cerr << "Found TITLE block." << endl;
                     given_title = true;
                     continue;     
                }
                if ( wordone == "element" ) {
                     isinblock = true;
                     inelement = true;
                     cerr << "Found ELEMENT block." << endl;
                     given_elements = true;
                     continue;     
                }               
                if ( wordone == "input" ) {
                     isinblock = true;
                     ininput = true;
                     cerr << "Found INPUT block." << endl;
                     given_input = true;
                     continue;     
                } 
                if ( wordone == "bond" ) {
                     isinblock = true;
                     inbond = true;
                     cerr << "Found BOND block." << endl;
                     given_bond = true;
                     bondwithin = defaultbondwithin; // in case it was set differently somewhere else!
                     //introduce WITHIN check here!
                     if ( unhashed.size() > 2 ) {
                          if ( unhashed.at(1) == "within" ) {
                               char dig = unhashed.at(2).at(0); // first char of third word
                               if ( isdigit( dig ) ){
                                    bondwithin = atof( unhashed.at(2).c_str() );
                                    cerr << "This block: bondwithin set to " << bondwithin << endl;
                               }
                          }
                     }
                     continue;     
                } 
                if ( wordone == "poly" ) {
                     isinblock = true;
                     inpoly = true;
                     cerr << "Found POLY block." << endl;
                     given_poly = true;
                     continue;     
                }
                if ( wordone == "option" ) {
                     isinblock = true;
                     inoption = true;
                     cerr << "Found OPTION block." << endl;
                     given_option = true;
                     continue;     
                }
                if ( wordone == "output" ) {
                     isinblock = true;
                     inoutput = true;
                     cerr << "Found OUTPUT block." << endl;
                     given_output = true;
                     continue;     
                }

           }
     
     
     } // done with WHILE on inputfile
  
  inputfile.close();
  return true;   
}

string tokenfromvar( int var ){
	if ( var == 0 ) return "0";
	if ( var == 1 ) return "+";
	if ( var == 2 ) return "-";
	return "0"; // trap case
}
int dirfromvar( int var ){
	if ( var == 0 ) return 0;
	if ( var == 1 ) return 1.;
	if ( var == 2 ) return -1.;
	return 0; // trap case	
}

string randomtag(){
	vector< string > tag;
	
	tag.push_back( "Don't type until you see the dots of their i's." );
	tag.push_back( "Program may fail in the event of Apocalypse." );
	tag.push_back( "Now wash your hands." );
	tag.push_back( "Please adjust your dress before leaving." );
	tag.push_back( "The first rule of Trappism is you do not talk about Trappism." );
	tag.push_back( "I'm afraid I can't let you do that, Dave." );
	tag.push_back( "Last day to send in your dollar!" );
	tag.push_back( "Splotch exists on artwork throughout." );
	tag.push_back( "I will keep your secret, but I require fine pastries for my silence." );
	tag.push_back( "Shall I eliminate him, Dance Captain?" );
	tag.push_back( "Ninja Professor. Silent. Deadly. Informative." );
	tag.push_back( "Neutron scattering is one of my hobbies." );
	tag.push_back( "And now, the end is near... and so I <ack>" );
	tag.push_back( "Your boss is standing behind you." );
	tag.push_back( "Empirical interatomic pistachio shell model." );
	tag.push_back( "Self-rolling dice render nerds obsolete." );
	tag.push_back( "Do not, I beg of you, install that which you cannot uninstall." );
	tag.push_back( "Consummatum est." );
	tag.push_back( "I am going to, or am about to, die - either expression is correct." );
	tag.push_back( "Who where what which why would you do that" );
	tag.push_back( "No fighting in the War Room, gentlemen." );
	tag.push_back( "Look what your careless hands have wrought." );
	tag.push_back( "Ia! Cthulhu f'thagn!" );
	
	int n = tag.size();
	int which = genrand_real2() * n; // 0 to n-1
	return tag.at(which );
}
