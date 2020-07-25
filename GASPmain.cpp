#include "GASPcontrol.h"


int main(int argc, char *argv[]){
    
  cerr << "Welcome to GASP version 6, for zeolites, MOFs and other frameworks," << endl;
  cerr << "   Shattered into a thousand shards, each with their own header file." << endl;
  cerr << "   Son of GASP 5, cousin of FRODA, and mayor of a small village on the coast." << endl;
  cerr << "   Handles structure input and output, bond finding, geometric analysis and relaxation." << endl;
  cerr << "NEW FEATURES filling a long-felt want:" << endl;
  cerr << "   Automation (partially) of flexibility window search;" << endl;
  cerr << "   CIF file support for input/output." << endl;
  cerr << "Written by Stephen Wells at Bath Uni." << endl;
  cerr << endl;

  init_genrand( static_cast< unsigned > (time(0) ) );
  
  //cerr << "Args: " << argc << endl;
  if ( argc == 1 ) {
     cerr << "No filename given; seeking default file gasp.inp" << endl;
     inpfilename = "gasp.inp"; 
  }
  else{
     stringstream whatin;
     whatin << argv[1];
     inpfilename = whatin.str();
     cerr << "Reading input from file " << inpfilename << endl;
  }
        
  //Step 1: parse input file
  bool inputokay = readcommands (inpfilename );
  if (inputokay){
     cerr << "Input okay, continuing." << endl;
  }
  else {
       cerr << "Input not found. Goodbye, cruel world." << endl;
       cerr << randomtag() << endl;
       exit(1);
  }
  
  
  //Step 2: validate input file
  
  bool valid = validation();
  if ( !valid ) {
     cerr << "Stopping execution due to invalid input logic." << endl;
     cerr << "Cretans are liars, said Epimenides the Cretan." << endl;
     cerr << randomtag() << endl;
     exit(1);     
  }
  
  //Step 3: obtain structure and bonding
  //read initial_structure vector< Atom > and initial_cellparam from input_structure_name
  bool readok;
  if ( useXTL ){
	  readok = inputstructure.readxtl( input_structure_name );
  }
  else{
	  readok = inputstructure.readcif( input_structure_name );
  }
  
  if ( !readok ) {
     cerr << "Stopping execution due to invalid structure input." << endl;
     cerr << "I require flawless crystals." << endl;
     cerr << randomtag() << endl;
     exit(1);     
  }
  if ( inputstructure.atom.size() == 0 ){
     cerr << "Problem - no atoms in input!" << endl;
     cerr << randomtag() << endl;
     exit(1);
  }
  else {
       cerr << "Read " << inputstructure.atom.size() << " atoms." << endl;
  }
  inputstructure.do14 = do_r14; //passing values for later use
  inputstructure.r14scale = r14times; 
  bool isinit = inputstructure.initialise( gridmin);
  //do cell vectors, inverse vectors, and cartesian positions, also fill grid
  if ( !isinit ){
	  cerr << "Sadly, we have choked on an unknown element. Please call an alchemist." << endl;
	  cerr << randomtag() << endl;
	  exit(1); // die
  }

  //if bonding file given, read bonding
  //give atoms their bond lists here.
  if ( given_bond_input ){
     bool bok = readbonds( input_bonding_name , bondline );
     if ( !bok ) {
        cerr << "Stopping execution due to invalid bonding input. Better bondage please." << endl;
        cerr << randomtag() << endl;
        exit(1);     
     }
     inputstructure.bondsintoatoms( bondline );
  }

  if ( !given_bond_input && given_poly ){
	inputstructure.polybondfind( polyspec, polypad );
  }

  if ( !given_bond_input && given_bond ){
	inputstructure.bondbondfind( bondspec );
  }

  //put bonding into bondline format ready for output and manipulation if needed
  if ( ! given_bond_input ){
	if ( !nodeloc ){
		vector< string > hybrids;
		hybrids.push_back("c");
		hybrids.push_back("n");
		inputstructure.hybridcheck( hybrids );
	} // replace later with deloc option
	inputstructure.formbondlines( bondline ); // bondlines to output if we need to
  }

  //okay, all bonding done
  //parse neighbours for poly specs
  if ( given_poly ){ 
	inputstructure.parsepoly( polyspec );
  }
  cerr << "Found " << inputstructure.poly.size() << " polyhedra." << endl;

  //now the non-poly bond clusters
  //note cluster array already contains null entry for each poly!
  if ( given_bond ){ // note given bond is required even if we have given bond input... tidy later?
		vector< Cluster > candidate = inputstructure.candidatecluster( bondspec ); // initial bonding clusters
		cerr << "Found " << candidate.size() << " non-poly clusters." << endl;          
		vector< Cluster > goodcandidate;
     
		// next loop will change when we go to "deloc" instead of "nodeloc"
		if ( nodeloc ){
			goodcandidate = candidate; // no change
		}
		else{
			goodcandidate = inputstructure.Garibaldi( candidate ); // reunification
		}

		//finally put into cluster and ghost!
		for ( int i = 0 ; i < goodcandidate.size() ; i++ ){
			inputstructure.cluster.push_back( goodcandidate.at(i) );
			inputstructure.ghost.push_back( goodcandidate.at(i) );
		}
		candidate.clear();
		goodcandidate.clear(); // just housekeeping.
  }
  
  cerr << "Found " << inputstructure.cluster.size() << " total clusters." << endl;

  if ( imperfect){
       inputstructure.prepimperfectclustersandghosts();
  }
  else{
       inputstructure.prepclustersandghosts();
  }    
  
  //polys now exist as poly, as cluster, and as ghost
  //nonpolys exist as cluster and as ghost
  // poly-ghosts have perfect geometry
  //label atoms by inclust! and isvertex
  inputstructure.labelvertices();
  //now do fourth neighbour things if called upon
  if ( do_r14 ){
       inputstructure.find14(); // find and report 1-4 cases
  }
  //initial structure has defined bonding, cluster and poly geometry.
  //initial_structure goes into structure now ready to move.
  //(i) read a new structure (positions) if ordered
  //(ii) apply a new cell if ordered

  structure = inputstructure; // :) use structure subsequently for e.g. relax

  if (given_new_structure){
	    bool readok;
	    if ( useXTL ){
			readok = newstructure.readxtl( new_structure_name );			
		}
		else{
			readok = newstructure.readcif( new_structure_name );
		}
		
        //validate or stop;
        if ( !readok ) {
			cerr << "Stopping execution due to invalid structure input." << endl;
            cerr << "The new structure is a bad structure." << endl;
            cerr << "I require flawless crystals." << endl;
            cerr << randomtag() << endl;
            exit(1);     
            }
         if ( newstructure.atom.size() == 0 ){
            cerr << "Problem - no atoms in input!" << endl;
            cerr << randomtag() << endl;
            exit(1);
         }
         else {
              cerr << "Read " << newstructure.atom.size() << " atoms." << endl;
         }
       
       //validate or stop;
       if ( newstructure.atom.size() != inputstructure.atom.size() ){
            cerr << "Problem with atom mismatch!" << endl;
            cerr << "New structure has " << newstructure.atom.size() << " atoms;" << endl;
            cerr << "Initial structure has " << inputstructure.atom.size() << " atoms." << endl;
            cerr << "I said no camels, that's five camels, can't you count?" << endl;
            cerr << randomtag() << endl;
            exit(1);
       }
		structure.newstructure( newstructure ); // updates positions and cell
  }
  
  if ( given_new_cell && !useGradual ){
       cerr << "Applying new cell parameters:" << endl;
       structure.newcell( newcellparam );
       cerr << structure.cell.printcell() << endl;
  }
  
  //Step 4: relax framework if requested
  
  if ( given_relax_command && !do_window && !useGradual ){
	  if ( do_jiggle ){
		  structure.jiggle( jigsize );
	  }
	  
	  int docheb=0;
	  if ( chebyshev ){
		  docheb = 1;
		  if ( use_auto_chebyshev ){
			  docheb = 2; 
		  }
	  }
	  bool goodrelax = false;
	  bool window = false;
	  goodrelax=gasprelax( structure, window, regrid_after, dampclash, docheb, cheb );
	  if (goodrelax) cerr << "GASP relaxation completed." << endl;
	  if (!goodrelax) cerr << "GASP relaxation not properly completed." << endl;
	  if ( window) cerr << "Relaxation IN WINDOW." << endl;	  
  }
  
  //"gradual" option occurs here
  if ( useGradual ){
	  //carry out nGradual+1 steps
	  //stepping proportionally from start to "new cell"
	  //step zero is relaxing the input structure so it goes to output 0
	  //so at each step:
	  //set the structure's cell appropriately between oldcellparam and newcellparam
	  //then carry out exactly the steps from relax itself
	  //and afterwards, write out appropriately
	  int docheb=0;
	  if ( chebyshev ){
		  docheb = 1;
		  if ( use_auto_chebyshev ){
			  docheb = 2; 
		  }
	  }	  
	  //cheb is a universal setting
	  vector< bool > waswindow; //keep track of windowing!
	  
	  vector< double > oldparam, nowparam;
	  oldparam = structure.cell.param;
	  nowparam = oldparam; //just to set size
	  for ( int i = 0 ; i < nGradual+1 ; i++ ){
		  for ( int j = 0 ; j < 6 ; j++ ){
			  double scale = 1.*i / ( nGradual ); //how far have we got?
			  nowparam.at(j) = scale * newcellparam.at(j) + (1.-scale)* oldparam.at(j);
		  }
		  //got cell
		  structure.newcell( nowparam );
		  cerr << "Cell at step " << i << endl;
		  cerr << structure.cell.printcell() << endl;
		  if ( do_jiggle ){
			  structure.jiggle( jigsize );
		  }		  
		   bool goodrelax = false;
		   bool window = false;
		   goodrelax=gasprelax( structure, window, regrid_after, dampclash, docheb, cheb );
		   if (goodrelax) cerr << "GASP relaxation completed." << endl;
		   if (!goodrelax) cerr << "GASP relaxation not properly completed." << endl;
		   if ( window) cerr << "Relaxation IN WINDOW." << endl;		  
		   waswindow.push_back( window );
		   //write out the structure to gradualname.rider

             stringstream namestream;
             namestream << gradualName << "." << setfill('0') << setw(2) << i;
             if ( useXTL ){
				 namestream << ".xtl";
			 } 
			 else{
				 namestream << ".cif";
			 }
             string outname = namestream.str();
             cerr << "Writing to file: " << outname << endl;
             stringstream titlestream;
             titlestream << "GRADUAL STEP STRUCTURE " << setfill('0') << setw(2) << i ;
             if (window ){
				 titlestream << " IN WINDOW" ;
			 }
			 else{
				 titlestream << " OUT OF WINDOW";
			 }
             string outtitle = titlestream.str();
             if ( useXTL ){
				 structure.writextl( outname, outtitle, do_label );
			 }
			 else{
				 structure.writecif( outname, outtitle );
			 }
		   
	  }
	  
	  //and do a final output summary
	  cerr << "Final summary after gradual action" << endl;
	  cerr << "STEP IN WINDOW?" << endl;
	  for ( int i = 0 ; i < nGradual + 1 ; i++ ){
		  cerr << i << " ";
		  if ( waswindow.at(i) ){
			  cerr << "YES" << endl;
		  }
		  else{
			  cerr << "NO" << endl;
		  }
	  }
	  
	  cerr << "Completed gradual process, stopping." << endl;
  }
  
  
  
  if (do_window){

     string winlogfile;
     winlogfile="WINDOWLOG.txt";
     if ( do_outputwindow ){
        stringstream newlog;
        newlog << windowbasename << ".windowlog.txt";
        winlogfile=newlog.str();
     }
     ofstream win;
     win.open( winlogfile.c_str(), ios::out | ios::app );
     
                 
	  cerr << "Hi, I'm the Window Search Routine!" << endl;
	  cerr << "You have asked for window type " << windowtype << endl;
	  win << "You have asked for window type " << windowtype << endl;
	  vector< Vector > sidestepN;
	  vector< Vector > anglestepN;
	  Vector ss;
	  Vector as;
	  if (windowtype == "cub"){
		  cerr << "Preparing CUBIC search type:" << endl;
		  cerr << "Direction key: [+++000]" << endl;
		  win << "Preparing CUBIC search type:" << endl;
		  win << "Direction key: [+++000]" << endl;
		  ss = Vector( 1., 1., 1.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [---000]" << endl;
		  win << "Direction key: [---000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }
	  else if (windowtype == "aaa"){
		  cerr << "Preparing AAA search type:" << endl;
		  cerr << "Direction key: [+00000]" << endl;
		  win << "Preparing CUBIC search type:" << endl;
		  win << "Direction key: [+00000]" << endl;
		  ss = Vector( 1., 0., 0.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [-00000]" << endl;
		  win << "Direction key: [-00000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }	  
	  else if (windowtype == "bbb"){
		  cerr << "Preparing BBB search type:" << endl;
		  cerr << "Direction key: [0+0000]" << endl;
		  win << "Preparing CUBIC search type:" << endl;
		  win << "Direction key: [0+0000]" << endl;
		  ss = Vector( 0., 1., 0.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [0-0000]" << endl;
		  win << "Direction key: [0-0000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }
	  else if (windowtype == "ccc"){
		  cerr << "Preparing CCC search type:" << endl;
		  cerr << "Direction key: [00+000]" << endl;
		  win << "Preparing CUBIC search type:" << endl;
		  win << "Direction key: [00+000]" << endl;
		  ss = Vector( 0., 0., 1.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00-000]" << endl;
		  win << "Direction key: [00-000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }
	  else if (windowtype == "alp"){
		  cerr << "Preparing ALPHA search type:" << endl;
		  cerr << "Direction key: [000+00]" << endl;
		  win << "Preparing CUBIC search type:" << endl;
		  win << "Direction key: [000+00]" << endl;
		  ss = Vector( 0., 0., 0.);
		  as = Vector( 1., 0., 0.);
		  as *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [000-00]" << endl;
		  win << "Direction key: [000-00]" << endl;
		  as *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }	  
	  else if (windowtype == "bet"){
		  cerr << "Preparing BETA search type:" << endl;
		  cerr << "Direction key: [0000+0]" << endl;
		  win << "Preparing CUBIC search type:" << endl;
		  win << "Direction key: [0000+0]" << endl;
		  ss = Vector( 0., 0., 0.);
		  as = Vector( 0., 1., 0.);
		  as *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [0000-0]" << endl;
		  win << "Direction key: [0000-0]" << endl;
		  as *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }
	  else if (windowtype == "gam"){
		  cerr << "Preparing GAMMA search type:" << endl;
		  cerr << "Direction key: [00000+]" << endl;
		  win << "Preparing CUBIC search type:" << endl;
		  win << "Direction key: [00000+]" << endl;
		  ss = Vector( 0., 0., 0.);
		  as = Vector( 0., 0., 1.);
		  as *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00000-]" << endl;
		  win << "Direction key: [00000-]" << endl;
		  as *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }
	  else if (windowtype == "iso"){
          //identify largest of a b c = L
          //Vector is a/L b/L c/L
          double myA = structure.cell.param.at(0);
          double myB = structure.cell.param.at(1);
          double myC = structure.cell.param.at(2);

          double myL = myA; // set to a
          if ( myB > myL ) myL = myB;
          if ( myC > myL ) myL = myC;
          
		  cerr << "Preparing ISOMETRIC search type:" << endl;
		  cerr << "Direction key: [+++000]" << endl;
		  win << "Preparing ISOMETRIC search type:" << endl;
		  win << "Direction key: [+++000]" << endl;
		  ss = Vector( myA/myL , myB/myL, myC/myL );
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [---000]" << endl;
		  win << "Direction key: [---000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );		  
	  }	  
	  else if (windowtype == "tet" || windowtype == "hex" ){
		  if ( windowtype == "tet" ){
               cerr << "Preparing TETRAGONAL search type." << endl;
               win << "Preparing TETRAGONAL search type." << endl;
          }
		  if ( windowtype == "hex" ){
               cerr << "Preparing HEXAGONAL search type." << endl;
               win << "Preparing HEXAGONAL search type." << endl;
          }		  
		  cerr << "Direction key: [++0000]" << endl;
		  win << "Direction key: [++0000]" << endl;
		  ss = Vector( 1., 1., 0.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [--0000]" << endl;
		  win << "Direction key: [--0000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00+000]" << endl;
		  win << "Direction key: [00+000]" << endl;
		  ss = Vector( 0., 0., 1.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00-000]" << endl;
		  win << "Direction key: [00-000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as ); 
		  
		  cerr << "Direction key: [+++000]" << endl;
		  win << "Direction key: [+++000]" << endl;
		  ss = Vector( 1., 1., 1.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [---000]" << endl;
		  win << "Direction key: [---000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [++-000]" << endl;
		  win << "Direction key: [++-000]" << endl;
		  ss = Vector( 1., 1., -1.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [--+000]" << endl;
		  win << "Direction key: [--+000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );

	  }
	  else if (windowtype == "six" ){
		  cerr << "Preparing SIX search type." << endl;		  
		  cerr << "Direction key: [+00000]" << endl;
		  win << "Preparing SIX search type." << endl;		  
		  win << "Direction key: [+00000]" << endl;
		  ss = Vector( 1., 0., 0.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [-00000]" << endl;	  
		  win << "Direction key: [-00000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [0+0000]" << endl;		  
		  win << "Direction key: [0+0000]" << endl;
		  ss = Vector( 0., 1., 0.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [0-0000]" << endl;	  
		  win << "Direction key: [0-0000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00+000]" << endl;		  
		  win << "Direction key: [00+000]" << endl;
		  ss = Vector( 0., 0., 1.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00-000]" << endl;	  
		  win << "Direction key: [00-000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [000+00]" << endl;		  
		  win << "Direction key: [000+00]" << endl;
		  ss = Vector( 0., 0., 0.);
		  as = Vector( 1., 0., 0.);
		  as *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [000-00]" << endl;	  
		  win << "Direction key: [000-00]" << endl;
		  as *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [0000+0]" << endl;		  
		  win << "Direction key: [0000+0]" << endl;
		  ss = Vector( 0., 0., 0.);
		  as = Vector( 0., 1., 0.);
		  as *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [0000-0]" << endl;	  
		  win << "Direction key: [0000-0]" << endl;
		  as *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00000+]" << endl;		  
		  win << "Direction key: [00000+]" << endl;
		  ss = Vector( 0., 0., 0.);
		  as = Vector( 0., 0., 1.);
		  as *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [00000-]" << endl;	  
		  win << "Direction key: [00000-]" << endl;
		  as *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
      }	  
	  else if (windowtype == "rho" ){
		  cerr << "Preparing RHOMBOHEDRAL search type." << endl;		  
		  cerr << "Direction key: [+++000]" << endl;
		  win << "Preparing RHOMBOHEDRAL search type." << endl;		  
		  win << "Direction key: [+++000]" << endl;
		  ss = Vector( 1., 1., 1.);
		  as = Vector( 0., 0., 0.);
		  ss *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [---000]" << endl;
		  win << "Direction key: [---000]" << endl;
		  ss *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [000+++]" << endl;
		  win << "Direction key: [000+++]" << endl;
		  ss = Vector( 0., 0., 0.);
		  as = Vector( 1., 1., 1.);
		  as *= windowstep; // setting step scale
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [000---]" << endl;
		  win << "Direction key: [000---]" << endl;
		  as *= -1.; // opposite direction!
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [++++++]" << endl;
		  win << "Direction key: [++++++]" << endl;
		  ss = Vector( 1., 1., 1.);
		  as = Vector( 1., 1., 1.);
		  ss *= windowstep; // setting step scale
		  as *= windowstep;
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [------]" << endl;
		  win << "Direction key: [------]" << endl;
		  ss *= -1.; // opposite direction!
		  as *= -1.;
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [+++---]" << endl;
		  win << "Direction key: [+++---]" << endl;
		  ss = Vector( 1., 1., 1.);
		  as = Vector( -1., -1., -1.);
		  ss *= windowstep; // setting step scale
		  as *= windowstep;
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as );
		  cerr << "Direction key: [---+++]" << endl;
		  win << "Direction key: [---+++]" << endl;
		  ss *= -1.; // opposite direction!
		  as *= -1.;
		  sidestepN.push_back( ss );
		  anglestepN.push_back( as ); 
	  }	  
	  else if (windowtype == "ort"){
		  cerr << "Preparing ORTHORHOMBIC search type:" << endl;
		  win << "Preparing ORTHORHOMBIC search type:" << endl;
		  
		  int avar; // counters 0 1 2 in the loops
		  int bvar;
		  int cvar;
		  string atok; // 0 + - tokens for reporting
		  string btok;
		  string ctok;
		  int adir; // -1 0 +1 vector components
		  int bdir;
		  int cdir;
		  
		  //for ( betavar = 0; betavar < 3; betavar++){
			  for ( cvar = 0; cvar < 3; cvar++){
				  for ( bvar = 0; bvar < 3; bvar++ ){
					  for ( avar = 0; avar < 3; avar++){
						  if ( avar == 0 && bvar == 0 && cvar == 0 ) continue; // pass all zeroes
						  //use special routines to convert from var to tok and dir
						  atok = tokenfromvar( avar );
						  btok = tokenfromvar( bvar );
						  ctok = tokenfromvar( cvar );
						  
						  adir = dirfromvar( avar );
						  bdir = dirfromvar( bvar );
						  cdir = dirfromvar( cvar );
						  
						  //construct and report string and vector
						  cerr << "Direction key: [";
						  cerr << atok << btok << ctok << "000]" << endl;
						  win << "Direction key: [";
						  win << atok << btok << ctok << "000]" << endl;
						  ss = Vector( adir, bdir, cdir);
						  as = Vector( 0., 0., 0.);
						  ss *= windowstep;
						  sidestepN.push_back( ss );
						  anglestepN.push_back( as );					  
						  
					  }
				  }
			  }
		  //}
		  		  
	  }
	  else if (windowtype == "mon"){
		  int avar; // counters 0 1 2 in the loops
		  int bvar;
		  int cvar;
		  int betavar;
		  string atok; // 0 + - tokens for reporting
		  string btok;
		  string ctok;
		  string betatok;
		  int adir; // -1 0 +1 vector components
		  int bdir;
		  int cdir;
		  int betadir;
		  
		  for ( betavar = 0; betavar < 3; betavar++){
			  for ( cvar = 0; cvar < 3; cvar++){
				  for ( bvar = 0; bvar < 3; bvar++ ){
					  for ( avar = 0; avar < 3; avar++){
						  if ( avar == 0 && bvar == 0 && cvar == 0 && betavar == 0 ) continue; // pass all zeroes
						  //use special routines to convert from var to tok and dir
						  atok = tokenfromvar( avar );
						  btok = tokenfromvar( bvar );
						  ctok = tokenfromvar( cvar );
						  betatok = tokenfromvar( betavar );
						  
						  adir = dirfromvar( avar );
						  bdir = dirfromvar( bvar );
						  cdir = dirfromvar( cvar );
						  betadir = dirfromvar( betavar );
						  
						  //construct and report string and vector
						  cerr << "Direction key: [";
						  cerr << atok << btok << ctok << "0" << betatok << "0]" << endl;
						  win << "Direction key: [";
						  win << atok << btok << ctok << "0" << betatok << "0]" << endl;
						  ss = Vector( adir, bdir, cdir);
						  as = Vector( 0., betadir, 0.);
						  ss *= windowstep;
						  as *= windowstep;
						  sidestepN.push_back( ss );
						  anglestepN.push_back( as );					  
						  
					  }
				  }
			  }
		  }
		  
		  
	  }
	  else if (windowtype == "tri"){
		  int avar; // counters 0 1 2 in the loops
		  int bvar;
		  int cvar;
		  int alphavar;
		  int betavar;
		  int gammavar;
		  string atok; // 0 + - tokens for reporting
		  string btok;
		  string ctok;
		  string alphatok;
		  string betatok;
		  string gammatok;
		  int adir; // -1 0 +1 vector components
		  int bdir;
		  int cdir;
		  int alphadir;
		  int betadir;
		  int gammadir;
		  
		  for ( gammavar = 0 ; gammavar < 3 ; gammavar++ ){
			  for ( betavar = 0; betavar < 3; betavar++){
				  for ( alphavar = 0 ; alphavar < 3; alphavar++){
					  for ( cvar = 0; cvar < 3; cvar++){
						  for ( bvar = 0; bvar < 3; bvar++ ){
							  for ( avar = 0; avar < 3; avar++){
								  if ( avar == 0 && bvar == 0 && cvar == 0 && alphavar == 0 && betavar == 0 && gammavar == 0 ) continue; // pass all zeroes
								  //use special routines to convert from var to tok and dir
								  atok = tokenfromvar( avar );
								  btok = tokenfromvar( bvar );
								  ctok = tokenfromvar( cvar );
								  alphatok = tokenfromvar( alphavar );
								  betatok = tokenfromvar( betavar );
								  gammatok = tokenfromvar( gammavar );
								  								  
								  adir = dirfromvar( avar );
								  bdir = dirfromvar( bvar );
								  cdir = dirfromvar( cvar );
								  alphadir = dirfromvar( alphavar );
								  betadir = dirfromvar( betavar );
								  gammadir = dirfromvar( gammavar );
								  								  
								  //construct and report string and vector
								  cerr << "Direction key: [";
								  cerr << atok << btok << ctok << alphatok << betatok << gammatok;
								  cerr << "]" << endl;
								  win << "Direction key: [";
								  win << atok << btok << ctok << alphatok << betatok << gammatok;
								  win << "]" << endl;
								  ss = Vector( adir, bdir, cdir);
								  as = Vector( alphadir, betadir, gammadir);
								  ss *= windowstep;
								  as *= windowstep;
								  sidestepN.push_back( ss );
								  anglestepN.push_back( as );					  
								  
							  }
						  }
					  }
				  }
			  }
		  }
		  
		  
	  }	  
	  
	  cerr << "Prepared " << sidestepN.size() << " search directions." << endl;
	  win << "Prepared " << sidestepN.size() << " search directions." << endl;	  
	  //cerr << "Oh." << endl;
	  //cerr << "Oh." << endl;
	  
	  Structure newpoint;
	  Structure topgood;
	  Vector topgoodside;
	  Vector topgoodangle;
	  Vector topbadside;
	  Vector topbadangle;
	  bool isgood;
	  bool newgood;
	  bool foundbad;
	  int direction;
	  int stepmul;
	  int maxmul = 4; // step limiter
	  Vector sidevec;
	  Vector anglevec;
	  Vector sidedelta;
	  Vector angledelta;
	  Vector sidestep;
	  Vector anglestep;

	  //begin by relaxing the initial structure! as we must fail if it doesn't relax
	  int docheb=0;
	  if ( chebyshev ){
		  docheb = 1;
		  if ( use_auto_chebyshev ){
			  docheb = 2; 
		  }
	  }
	  bool goodrelax = false;
	  bool window = false;
	  goodrelax=gasprelax( structure, window, regrid_after, dampclash, docheb, cheb );
	  if (goodrelax){
             cerr << "GASP relaxation of starting structure completed." << endl;
             win << "GASP relaxation of starting structure completed." << endl;
      }
	  if (!goodrelax) cerr << "GASP relaxation of starting structure not properly completed." << endl;
	  if (!goodrelax) win << "GASP relaxation of starting structure not properly completed." << endl;
	  if ( window){
		   cerr << "Relaxation IN WINDOW with ";	
		   cerr << "starting cell: " << structure.cell.printcell() << " Volume: " << structure.cell.volume() << endl;
		   win << "Relaxation IN WINDOW with ";	
		   win << "starting cell: " << structure.cell.printcell() << " Volume: " << structure.cell.volume() << endl;
	  }
	  
      if ( !goodrelax || !window ){
		  cerr << "Cannot proceed with window search: starting structure not in window." << endl;
		  win << "Cannot proceed with window search: starting structure not in window." << endl;
		  cerr << randomtag() << endl;
		  exit(1); // I die now.
	  }
	  
	  for ( int s = 0 ; s < sidestepN.size(); s++ ){
		  //setup phase
		  newpoint = structure;
		  topgood = structure;
		  sidestep = sidestepN.at(s);
		  anglestep = anglestepN.at(s);
		  sidevec = first3( structure.cell.param );
		  anglevec = second3( structure.cell.param );
		  topgoodside = sidevec;
		  topgoodangle = anglevec;
		  topbadside = nullvec;
		  topbadangle = nullvec;
		  isgood = true; // initially relaxed if we've got here.
		  direction = 1;
		  stepmul = 1;
		  foundbad = false; // initially not found a bad one
		  bool isdone = false; // have we finished?
		  
		  cerr << "Preparing to search direction: " << sidestep << " " << anglestep << endl;
		  win << "Preparing to search direction: " << sidestep << " " << anglestep << endl;
          		  
		  do {
			  if ( foundbad ){
				  //stop trap first
				  sidedelta = topgoodside - topbadside;
				  angledelta = topgoodangle - topbadangle;
				  if ( ( (sidedelta.sq()- sidestep.sq()) < 1e-8 ) && ( (angledelta.sq()-anglestep.sq()) < 1e-8 ) ){
                       //the difference is one step, ergo we have found an edge
                       isdone = true; // time to stop!
                       double vol = topgood.cell.volume(); 
                       cerr << "FOUND WINDOW EDGE: " << topgoodside << " " << topgoodangle << " Volume: " << vol << endl;
                       win << "FOUND WINDOW EDGE: " << topgoodside << " " << topgoodangle << " Volume: " << vol << endl;
                       continue; // cycle loop to halt condition
                  }
				  
				  //proceed with binomial search
				  //update the values;
				  Vector change = sidestep * stepmul * direction;
				  sidevec += change;
				  change = anglestep*stepmul*direction;
				  anglevec += change;
				  cerr << "TRYING: " << sidevec << " " << anglevec << endl;
				  win << "TRYING: " << sidevec << " " << anglevec << "(" <<stepmul << "," << direction << ")" << endl;
				  vector< double > newpars = sixfromtwovecs( sidevec, anglevec );
				  newpoint.newcell( newpars ); // new cell params in newpoint structure
				  newgood = false;
				  bool goodrel = gasprelax( newpoint, newgood, regrid_after, dampclash, docheb, cheb );
				  if (!goodrel || !newgood ){
                     //relaxation has failed
                     newgood = false; // just to make sure
                  }
                  if ( stepmul > 1) stepmul /= 2; // binomial division of step size
                  if ( newgood != isgood ) direction *= -1; //change direction if we crossed the edge;
                  if (newgood ){
                     //this is a good relaxed point
                     topgood = newpoint; // best so far by construction
                     isgood = true;
                     double vol = topgood.cell.volume();
                     cerr << "GOOD point: " << sidevec << " " << anglevec << " Volume: " << vol << endl;
                     win << "GOOD point: " << sidevec << " " << anglevec << " Volume: " << vol << endl;
                     topgoodside = sidevec;
                     topgoodangle = anglevec;
                  }
				  else{
                       //this is a bad point
                       //recover from bad if need be:
                       if (!goodrel) newpoint = topgood; // we''ll update cell anyway
                       isgood = false;
                       cerr << "BAD point: " << sidevec << " " << anglevec << endl;
                       win << "BAD point: " << sidevec << " " << anglevec << endl;
                       topbadside = sidevec;
                       topbadangle = anglevec;
                  }
				  
			  }
			  else{
				  //not found bad yet
				  //update the values;
				  Vector change = sidestep * stepmul * direction;
				  sidevec += change;
				  change = anglestep*stepmul*direction;
				  anglevec += change;				  
				  cerr << "TRYING: " << sidevec << " " << anglevec << endl;
				  win << "TRYING: " << sidevec << " " << anglevec << "(" <<stepmul << "," << direction << ")" << endl;
				  vector< double > newpars = sixfromtwovecs( sidevec, anglevec );
				  newpoint.newcell( newpars ); // new cell params in newpoint structure
				  newgood = false; // checking
				  bool goodrel = gasprelax( newpoint, newgood, regrid_after, dampclash, docheb, cheb );				  
				  if ( !goodrel || !newgood ){
                     //relaxation has failed
                     newgood = false; // just to make sure
                     foundbad = true;
                  }				  
                  if ( newgood ){
                     //this is a good relaxed point
                     topgood = newpoint;
                     isgood = true;
                     double vol = topgood.cell.volume();          
                     cerr << "GOOD point: " << sidevec << " " << anglevec << " Volume: " << vol << endl;
                     win << "GOOD point: " << sidevec << " " << anglevec << " Volume: " << vol << endl;
                     topgoodside = sidevec;
                     topgoodangle = anglevec;
                     if ( stepmul < maxmul) stepmul *= 2; // take bigger steps
                  }
				  else{
                       //this is a bad point
                       //recover from bad if need be:
                       if (!goodrel) newpoint = topgood; // we''ll update cell anyway
                       isgood = false;
                       cerr << "BAD point: " << sidevec << " " << anglevec << endl;
                       win << "BAD point: " << sidevec << " " << anglevec << endl;
                       topbadside = sidevec;
                       topbadangle = anglevec;
                       if ( stepmul > 1 ) stepmul /= 2; // reduce step...
                       direction = -1; // and reverse direction
                  }                  
				  
			  }	  
			  
		  } while (!isdone ); // loop ends
		  //we can report "topgood" structure now if desired
		  
		  cerr << "Finished search direction: " << sidestep << " " << anglestep << endl;
          cerr << "Finished search " << s+1 << " of " << sidestepN.size() << endl;		  
		  win << "Finished search direction: " << sidestep << " " << anglestep << endl;
          win << "Finished search " << s+1 << " of " << sidestepN.size() << endl;
          if ( do_outputwindow ){
             stringstream edgenamestream;
             edgenamestream << windowbasename << "." << setfill('0') << setw(2) << s+1;
             if ( useXTL ){
				 edgenamestream << ".xtl";
			 } 
			 else{
				 edgenamestream << ".cif";
			 }
             string edgename = edgenamestream.str();
             win << "Writing to file: " << edgename << endl;
             stringstream edgetitlestream;
             edgetitlestream << "WINDOW EDGE STRUCTURE " << setfill('0') << setw(2) << s+1 ;
             string edgetitle = edgetitlestream.str();
             if ( useXTL ){
				 topgood.writextl( edgename, edgetitle, do_label );
			 }
			 else{
				 topgood.writecif( edgename, edgetitle );
			 }
             
          }

	  }
	  
  }


  //Step 5: file output
  if ( do_output_structure && useXTL ){
	  structure.writextl( output_structure_name, title, do_label );
  }
  else if ( do_output_structure ){
	  structure.writecif( output_structure_name, title );
  }
  if ( do_output_bonding ){
       writebonding( output_bonding_name, bondline ); // method is in Geometry
  }
  if ( do_output_poly ){
       structure.writepol( output_poly_name );
  }  
  if ( do_output_angles ){
       structure.writeangles( output_angles_name );
  }  
  if ( do_output_bonds ){
       structure.writebondlengths( output_bonds_name );
  }
  if ( do_output_mismatches ){
       structure.writemismatches(output_mismatches_name );      
  }

    cerr << randomtag() << endl;
    exit(0); //final exit
} // END OF MAIN


