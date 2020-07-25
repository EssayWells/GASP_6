#include "Structure.h"

string Cell::printcell (){
       stringstream reportcell;
       //reportcell << "Cell a-b-c-alpha-beta-gamma:";
       for ( int i = 0; i < 6 ; i++ ) {
           reportcell << " " << fixed << showpoint << setw(8) << setprecision(4) << param.at(i);
       }
       return reportcell.str();
}

double Cell:: volume(){
       Vector foo;
       foo= cellj.cross( cellk );
       double vol = celli.dot( foo );
       return vol;
}

//
//bonded returns true if atoms a1 and a2 have any clust in common
//
bool bonded( Atom &a1, Atom &a2){
     for ( int i = 0 ; i < a1.inclust.size() ; i++ ){
         for ( int j = 0 ; j < a2.inclust.size() ; j++ ){
             if ( a1.inclust.at(i) == a2.inclust.at(j) ) return true; // common cluster, bonded
         }
     }
     return false;
}




//
//makecell creates the celli cellj cellk vectors from the cell parameters
//
void Cell::makecell(){
  cerr << "Making cell vectors from parameters." << endl;
  double twopi360 = mytwopi/360.0;
  double radalpha, radbeta, radgamma;

  radalpha = twopi360 * param.at(3);
  radbeta = twopi360 * param.at(4);
  radgamma = twopi360 * param.at(5);

  cellj = Vector( 0.0, param.at(1), 0.0); // aligned on y axis
  cellk = Vector( 0.0, param.at(2)*cos(radalpha), param.at(2)*sin(radalpha) ); // lies in yz plane
  celli.y = param.at(0)*cos(radgamma);
  celli.z = param.at(0)*( cos(radbeta) - cos(radgamma) * cos(radalpha) ) / sin(radalpha);
  celli.x = sqrt( pow( param.at(0),2) - pow( celli.y,2) -pow(celli.z,2) ); 
  cerr << "Cell vectors made." << endl;
  cerr << celli << endl;
  cerr << cellj << endl;
  cerr << cellk << endl;
  return;
}

//
//makestar creates reciprocal vectors matching the cell vectors.
//
void Cell::makestar(){
  cerr << "Making reciprocal vectors." << endl;
  Vector dummyvec;
  double dot;

  dummyvec = celli.cross(cellj);
  dot = dummyvec.dot( cellk );
  kstar = dummyvec/dot;
  
  dummyvec =  cellj.cross(cellk);
  dot = dummyvec.dot(  celli);
  istar = dummyvec/dot;

  dummyvec = cellk.cross(celli);
  dot = dummyvec.dot(  cellj);
  jstar = dummyvec/dot;
  // star vectors dot to delta with cell vectors

  cerr << "Reciprocal vectors made." << endl;
  cerr << istar << endl;
  cerr << jstar << endl;
  cerr << kstar << endl;
  return;
} //ENDOFMAKESTAR

void Cell::makegrid( ){ // creates coarse grid
  
  if (mygridmin < 3.0 ) mygridmin = 3.0; // just a trap, lose if not needed
  // grid sizing logic
  double aspace = 1.0/(sqrt(  istar.sq() ) );
  double bspace = 1.0/(sqrt(  jstar.sq() ) ); 
  double cspace = 1.0/(sqrt(  kstar.sq() ) );   
  Nx = (int) ( aspace / mygridmin );
  Ny = (int) ( bspace / mygridmin );
  Nz = (int) ( cspace / mygridmin );

  grid.clear();
  grid.resize(Nx*Ny*Nz);
  
  cerr << "Coarse gridding: " << Nx << "," << Ny << "," << Nz << " boxes." << endl;

  return;
}


//
// fractocart uses cell to turn frac position into Cartesian position
//
Vector Cell::fractocart( Vector& fracin ){
       Vector result;
       result = celli * fracin.x + cellj * fracin.y + cellk * fracin.z;
       return result;
}

//
//carttofrac uses *star to turn Cartesian position into fractional coordinate
//
Vector Cell::carttofrac( Vector& cartin ){
       Vector result;
       result.x = cartin.dot( istar );
       result.y = cartin.dot( jstar );
       result.z = cartin.dot( kstar );
       return result;
}
       
//
//wrapcartpos ensures that a Cartesian position lies within the first positive unit cell.
//       
Vector Cell::wrapcartpos( Vector& cartin ){
       Vector resultA;
       resultA = carttofrac ( cartin ) ;//now it's a fracpos
       Vector resultB;
       resultB = wrapfracpos( resultA) ;// wrapped fractional
       resultA = fractocart( resultB) ; // and now it's a cart again
       return resultA;
}

//
//wrapcartdel ensures that a Cartesian delta lies within half a unit cell of the origin
//
Vector Cell::wrapcartdel( Vector& cartin ){
       Vector resultA;
       resultA = carttofrac ( cartin ) ;//now it's a fracdel
       Vector resultB;
       resultB = wrapfracdel( resultA) ;// wrapped fractional
       resultA = fractocart( resultB) ; // and now it's a cart again
       return resultA;
}


//
//routines to return grid info for a given position
//
int Cell::posNx( Vector &pos ){
	double f = pos.dot( istar ); // make frac x
	int result = f * Nx;
	if ( result == Nx ) result--;
	return result;
}
int Cell::posNy( Vector &pos ){
	double f = pos.dot( jstar ); // make frac y
	int result = f * Ny;
	if ( result == Ny ) result--;
	return result;	
}
int Cell::posNz( Vector &pos ){
	double f = pos.dot( kstar ); // make frac z
	int result = f * Nz;
	if ( result == Nz ) result--;
	return result;
}
int Cell::fracNx( double fx ){
	int result = fx * Nx;
	if ( result == Nx ) result--;
	return result;
}
int Cell::fracNy( double fy ){
	int result = fy * Ny;
	if ( result == Ny ) result--;
	return result;	
}
int Cell::fracNz( double fz ){
	int result = fz * Nz;
	if ( result == Nz ) result--;
	return result;	
}
int Cell::Ngrid( int ax, int ay, int az){
	int result = az + ay*Nz + ax*Ny*Nz;
	return result; // the "whichgrid" to look in
}

void Cell::remove( int index, int gx, int gy, int gz ){
	//cerr << "In cell.remove" << endl;
	int where = Ngrid( gx, gy, gz ); // look in this cell
	//cerr << "Looking in cell " << where << endl;
	for ( int q = grid.at( where ).size() - 1; q > -1 ; q-- ){
		//cerr << "Trying entry " << q << endl;
		if ( grid.at( where ).at(q) == index ){
			grid.at(where).erase( grid.at(where).begin() +q ); //kill the right item
			//cerr << "Found entry at index " << q << endl;
			break; // and we're done
		}
	}
}
void Cell::replace( int index, int gx, int gy, int gz ){
	int where = Ngrid( gx, gy, gz ); // look in this cell
	grid.at( where ).push_back( index ); // I hope you mean it
}

// take atom index aid out of cell grid
void Structure::gridremove( int aid ){
	//cerr << "In gridremove." << endl;
	cell.remove( aid, atom.at(aid).gridx, atom.at(aid).gridy, atom.at(aid).gridz );
	//cerr << "Leaving gridremove." << endl;
} 
// put index aid into cell grid
void Structure::gridreplace( int aid ){
	//cerr << "In gridreplace." << endl;
	cell.replace( aid, atom.at(aid).gridx, atom.at(aid).gridy, atom.at(aid).gridz );
	//cerr << "Leaving gridremove." << endl;
}

//
//fillgrid places all atoms into cells in the coarse grid.
//
//only run after makegrid has produced a clean empty grid
bool Structure::fillgrid ( ){
	 
  int gridX,gridY,gridZ;
  int whichgrid;

  for (int i =0; i < atom.size(); i++){
	  
    gridX = (int) (atom.at(i).fracpos.x * cell.Nx);
    gridY = (int) (atom.at(i).fracpos.y * cell.Ny);
    gridZ = (int) (atom.at(i).fracpos.z * cell.Nz);

    if ( gridX < 0 || gridY < 0 || gridZ < 0 ){
         cerr << "ERROR: Atom " << i+1 << " grid " << gridX << " " << gridY << " " << gridZ << endl;
         cerr << "Atom " << i+1 << " is at " << atom.at(i).fracpos.x << "," << atom.at(i).fracpos.y << "," << atom.at(i).fracpos.z << endl;
         return false; // flags grid failure
    }
    if (gridX >= cell.Nx ) gridX = cell.Nx -1 ; // this traps numerical slop at the very top edge of the cell.
    if (gridY >= cell.Ny ) gridY = cell.Nx -1 ;
    if (gridZ >= cell.Nz ) gridZ = cell.Nx -1 ;

    whichgrid = gridZ + gridY*cell.Nz + gridX*cell.Ny*cell.Nz;
    if ( whichgrid < 0 || whichgrid >=cell.grid.size() ){
         cerr << "ERROR: Atom " << i+1 << " whichgrid " << whichgrid << endl;
         cerr << "Atom " << i+1 << " is at " << atom.at(i).fracpos.x << "," << atom.at(i).fracpos.y << "," << atom.at(i).fracpos.z << endl;         
         return false; // flags grid failure
    }
    //cerr << "Placing atom " << i << " in grid " << gridX << "," << gridY << "," << gridZ << ", so " << whichgrid << endl;
    cell.grid.at(whichgrid).push_back(i);
    atom.at(i).gridx = gridX;
    atom.at(i).gridy = gridY;
    atom.at(i).gridz = gridZ;
  }
  return true; // good grid fill
  //cerr << "Grid filled." << endl;
}


//
//formcluster sets up the cluster cpos and bonds
//
void Structure::formcluster( Cluster &clust ){
     vector< Vector > position;
     Vector tempvec;
     position.clear();
     Vector workvec;
     
     for ( int i = 0; i < clust.members.size(); i++ ){
         position.push_back( atom.at( clust.members.at(i) ).pos );
     }
     //start with anchor atom.
     workvec = position.at(0);
     for ( int i = 1; i < clust.members.size(); i++ ){
         tempvec = position.at(i) - workvec; // bond vector, raw
         position.at(i) = cell.wrapcartdel( tempvec ) ;// bond vector, wrapped
         position.at(i) += workvec; // new cartesian including wrapping
     }
     //now we have a floating set of wrapped Cartesian positions
     workvec = Vector(0,0,0); // nulling
     for ( int i = 0; i < clust.members.size(); i++ ){
         workvec += position.at(i); // sum of positions
     }
     workvec /= clust.members.size(); // average!
     clust.cpos = workvec; // centre position done!
     //cerr << "Cluster centre: " << clust.cpos << endl;
     for ( int i = 0; i < clust.members.size(); i++ ){
         position.at(i) -= clust.cpos ; // now they're relative vectors!
         //double d2 = sqrt( position.at(i).sq() ) ;
         //cerr << "Cluster vector: " << position.at(i) << " , (" << d2 << ")" << endl;
     }
     clust.bond = position; // done?     
     return;
}

//
//
//
void Structure::writextl( string filename, string title, bool do_label /* = false */ ){
     cerr << "Writing structure to " << filename << endl;
     ofstream oxtlfile( filename.c_str() );

     oxtlfile << "TITLE " << title << endl;

     oxtlfile << "CELL" << endl;
     oxtlfile << cell.printcell() << endl;
     oxtlfile << "ATOMS " << endl;
     oxtlfile << "NAME X Y Z" << endl;
     
     for ( int i = 0 ; i < atom.size(); i++ ){
         oxtlfile << atom.at(i).species << " ";
         oxtlfile << atom.at(i).fracpos;
         if ( do_label ){
              oxtlfile << " " << i+1;
         }
         oxtlfile << endl;
     }

     
     oxtlfile << "EOF" << endl;
     oxtlfile.close();   
     cerr << "Done writing to file " << filename << endl;       
}


//
//readxtl pulls xtl data out of the given file and into the structure arrays
//
bool Structure::readxtl( string filename ){
     //read initial_structure vector< Atom > and initial_cellparam from input_structure_name
     cerr << "Seeking to read structure from " << filename << endl;
     atom.clear();
     cell.param.clear();
     cell.param.resize(6);
     
     Vector fracpos;
     Atom dummyatom;
     bool reading_atoms = false;
     bool reading_cell = false;
     bool haseof = false;
     bool hascell = false;
     bool hasatoms = false;
     
     string linein;
     string foo; //use as buffer
     vector< string > tokensin;
     
     ifstream xtlfile( filename.c_str() );
     if ( xtlfile.fail() ) {
        cerr << "Problem! Cannot read file " << filename << endl;
        return false;     
     }

     while( getline(xtlfile, foo) ){
            dummyatom = nullatom;
            linein = to_lower(foo);// lose case for simplicity
            //cerr << linein << endl;
            tokensin = chop(linein, " \t");
            
            if (tokensin.size() == 0 ) continue; //blank line?
            //for ( int k = 0; k < tokensin.size() ; k++ ){
			//	cerr << k << "/" << tokensin.at(k) << ";";
			//}
			//cerr << endl;
            
            
            if (tokensin.at(0) == "eof" ) {
               haseof = true; // we read an eof line
               cerr << "EOF" << endl;
               break; // reached end of xtl format
            }
            else if ( tokensin.at(0) == "cell" ) {
                 reading_cell = true;
                 cerr << "CELL" << endl;
                 continue; // next line should be cell params
            }
            else if ( tokensin.at(0) == "name" ) {
                 reading_atoms = true;
                 hasatoms = true;
                 cerr << "ATOMS" << endl;
                 continue; //
            }            
            else if ( reading_atoms ) {
                 if ( tokensin.size() < 4 ) {
                      cerr << "Problem reading from line " << linein << endl;
                      cerr << "Bad atom line, needs element and frac-x y z" << endl;
                      return false; // choked
                 }
                 dummyatom = nullatom;
                 dummyatom.species = tokensin.at(0);
				 // element check will be done in initialise
                 fracpos.x = atof(tokensin.at(1).c_str());  
                 fracpos.y = atof(tokensin.at(2).c_str());  
                 fracpos.z = atof(tokensin.at(3).c_str());  
                 dummyatom.fracpos = wrapfracpos(fracpos); // read into unit cell range
                 dummyatom.role = 0;
                 dummyatom.mobile = true;
                 dummyatom.hybrid = 0 ;// generic case; check for sp2 later
                 atom.push_back(dummyatom);
            }
            else if ( reading_cell ) {
                 if ( tokensin.size() < 6 ) {
                      cerr << "Problem reading from line " << linein << endl;
                      cerr << "Bad cell line, needs a b c alpha beta gamma" << endl;
                      return false; // choked
                 }
                 //cerr << endl;
                 for (int i=0; i<6; i++){
					 //cerr << tokensin.at(i) << " ";
                     cell.param.at(i) = atof(tokensin.at(i).c_str());
                     //cerr << cell.param.at(i) << endl;
                 }
                 reading_cell = false; // otherwise it will choke on the next line
                 hascell=true;
                 cerr << cell.printcell() << endl;
            }
            
     } // end of while
    
     if ( haseof && hascell && hasatoms ){
          return true; // all good
     }
     return false;
}


//
//nearlist returns the IDs of all atoms near to atom "me" according to coarse grid
//
vector < int > Structure::nearlist( int me, double within ){
       vector< int > result;
       result.clear();
       
       //grid IDs;
       int tg, og ; // this grid cell, other grid cell
       int tx,ty,tz; // grid cell x,y,z occupied by "me"
       int ox,oy,oz; // grid cell to check
       
       tx = atom.at(me).gridx;
       ty = atom.at(me).gridy;
       tz = atom.at(me).gridz;
       
       Vector diff; // initial difference before wrapping
       Vector deltaf;// fractional difference in position
       Vector delta; // Cartesian difference
       double d2; // distance-squared
       double threshold = within * within; // use distance-squared to avoid sqrt cost
       
       for ( int dx = -1 ; dx < 2 ; dx++ ){
           ox = (tx + dx );
           if ( ox == -1 ) ox += cell.Nx; // mind negative!
           if ( ox == cell.Nx ) ox = 0 ; // top of cell
           for ( int dy = -1; dy < 2 ; dy++ ){
               oy = ( ty + dy ) % cell.Ny ;
               if ( oy == -1 ) oy += cell.Ny; // mind negative!
               if ( oy == cell.Ny ) oy = 0 ; // top of cell               
               for ( int dz = -1 ; dz < 2 ; dz++ ){
                   oz = ( tz + dz ) % cell.Nz;
                   if ( oz == -1 ) oz += cell.Nz; // mind negative!
                   if ( oz == cell.Nz ) oz = 0 ; // top of cell                   
                   og = oz + ( oy * cell.Nz ) + ( ox * cell.Ny * cell.Nz );
                   //cerr << "Looking in grid " << og << " of " << grid.size() << endl;
                   // look in og for possible neighbours
                   for ( int j = 0; j < cell.grid.at(og).size(); j++){
                       int you = cell.grid.at(og).at(j); // other atom!
                       if ( me == you ) continue; // Malkovich
                       diff = atom.at(me).fracpos - atom.at(you).fracpos;
                       deltaf = wrapfracdel ( diff ) ; // shortest vector, not across cell
                       delta = cell.fractocart( deltaf) ; // Cartesians now
                       d2 = delta.dot( delta ); // distance-squared!
                       if ( d2 < threshold ){
                            //unexpected bug: in a small cell this can find the same atom twice! once wrapped and once not
                            bool already = false;
                            for ( int k = 0 ; k < result.size(); k++) {
                                if ( you == result.at(k) ){
                                   already = true; // got a copy already
                                   break; // we know to leave now
                                }
                            }
                            if ( already ) continue; // go to next round of loop
                            result.push_back( you );
                            //cerr << "Atom " << me << " near atom " << you << endl;
                       }
                   }
                  
               }
           }
       }
       return result; // vector of neighbour ids!
}



//
//fullnearlist returns the nearlist for all atoms in the structure, avoiding double count
//
vector< vector < int > > Structure::fullnearlist( double within ){
        vector< vector < int > > result;
        result.clear();
        
        int nat = atom.size();
        
        result.resize( nat ); // one array per atom
        
        for ( int i = 0 ; i < nat ; i++ ){
            int tg, og ; // this grid cell, other grid cell
            int tx,ty,tz; // grid cell x,y,z occupied by "me"
            int ox,oy,oz; // grid cell to check
       
            tx = atom.at(i).gridx;
            ty = atom.at(i).gridy;
            tz = atom.at(i).gridz;
       
            Vector diff; // initial difference before wrapping
            Vector deltaf;// fractional difference in position
            Vector delta; // Cartesian difference
            double d2; // distance-squared
            double threshold = within * within; // use distance-squared to avoid sqrt cost
       
            for ( int dx = -1 ; dx < 2 ; dx++ ){
                ox = (tx + dx );
                if ( ox == -1 ) ox += cell.Nx; // mind negative!
                if ( ox == cell.Nx ) ox = 0 ; // top of cell
                for ( int dy = -1; dy < 2 ; dy++ ){
                    oy = ( ty + dy ) % cell.Ny ;
                    if ( oy == -1 ) oy += cell.Ny; // mind negative!
                    if ( oy == cell.Ny ) oy = 0 ; // top of cell               
                    for ( int dz = -1 ; dz < 2 ; dz++ ){
                        oz = ( tz + dz ) % cell.Nz;
                        if ( oz == -1 ) oz += cell.Nz; // mind negative!
                        if ( oz == cell.Nz ) oz = 0 ; // top of cell                   
                        og = oz + ( oy * cell.Nz ) + ( ox * cell.Ny * cell.Nz );
                        //cerr << "Looking in grid " << og << " of " << grid.size() << endl;
                        // look in og for possible neighbours
                        for ( int j = 0; j < cell.grid.at(og).size(); j++){
                            int you = cell.grid.at(og).at(j); // other atom!
                            if ( i >= you ) continue; // Malkovich
                            diff = atom.at(i).fracpos - atom.at(you).fracpos;
                            deltaf = wrapfracdel ( diff ) ; // shortest vector, not across cell
                            delta = cell.fractocart( deltaf) ; // Cartesians now
                            d2 = delta.sq() ; // distance-squared!
                            if ( d2 < threshold ){
								//unexpected bug: in a small cell this can find the same atom twice! once wrapped and once not
								bool already = false;
								for ( int k = 0 ; k < result.at(i).size(); k++) {
									if ( you == result.at(i).at(k) ){
										already = true; // got a copy already
										break; // we know to leave now
									}
								}
								if ( already ) continue; // go to next round of loop								
								
                               result.at( i ).push_back( you );

                            }
                        }
                  
                  }
                }
            }            
            
            
            
       }
        
        return result;       
}


//
//mymis returns the vector from Atom to its vertex in Cluster
//
Vector Structure::mymis( Atom &at , Cluster &clus, int vert ){
       Vector result;
       Vector diff;
       //atom, cluster, which vertex of the cluster
       diff = clus.cpos + clus.bond.at(vert); // Cartesian vertex position;
       diff -= at.pos; // Cartesian difference from at to vert...
       result = cell.wrapcartdel ( diff ); // allow for cell wrapping!
       
       return result;
}

//
//hybridcheck checks species and bonding to establish hybrid status
//
void Structure::hybridcheck( vector< string > hybridel ){
	for ( int i =0; i < atom.size(); i++) {
		if ( atom.at(i).bondto.size() > 3 ) continue; // not sp3
		for ( int j = 0; j < hybridel.size(); j++) {
			string atsp = to_lower( atom.at(i).species );
			string hyb = to_lower( hybridel.at(j) );
			if ( atsp == hyb ){
					atom.at(i).hybrid = j+1; // is the j+1th hybrid type
			}
		}
	}	
}

void Structure::mcrigidcheck(){
	for ( int i = 0 ; i < atom.size(); i++ ){
		//initial label step: all hybrid
		atom.at(i).hybrid = 1;
	}
	//now remove label from metal centres and their neighbours
	for ( int i = 0 ; i < atom.size(); i++ ){
		int whichr = atom.at(i).role;
		if ( whichr != -1 ) continue; // skip non-metal-centres
		atom.at(i).hybrid = 0;
		for ( int j = 0; j < atom.at(i).bondto.size(); j++ ){
			int other = atom.at(i).bondto.at(j);
			atom.at(other).hybrid = 0; // unhybrid vertices
		}
	}
}

bool Structure::elecheck( Atom &atom ){
		atom.type = -1;
		for ( int e = 0 ; e < element.size() ; e++ ){
			if ( lowerCaseVersion( atom.species ) == lowerCaseVersion( element.at(e).species ) ){
				atom.type = e;
                atom.radius = element.at(e).radius;
                atom.role = element.at(e).role; // transfers role for MC use
                break;
            }
        }
        if ( atom.type == -1 ){
			cerr << "Unrecognised element: " << atom.species << endl;
			return false; // choked
		}	
		return true;
}

int Structure::eletype( string &sp ){
	int result = -1;
	for ( int e = 0 ; e < element.size(); e++ ){
		if ( lowerCaseVersion( sp ) == lowerCaseVersion( element.at(e).species ) ){
			result = e;
			break;
		}
	}
	return result;
}

bool Structure::initialise( double thisgridmin ){
	cell.makecell();
	cell.makestar();
	
	for ( int i =0; i < atom.size(); i++) {
		atom.at(i).pos = cell.fractocart( atom.at(i).fracpos );
		atom.at(i).initial_pos = atom.at(i).pos;
		// do element check here as well:
		bool gotel = false;
		gotel = elecheck( atom.at(i) );
		if ( !gotel ) return false; //choked on bad element	
	}
	cerr << "Scanned " << atom.size() << " atoms for position and role." << endl;
	cell.mygridmin = thisgridmin;
	cell.makegrid();
	fillgrid();
	return true;
}

void Structure::prepclustersandghosts(){
	
	for ( int i = 0; i < cluster.size(); i++ ){
      formcluster( cluster.at(i) ); // checking
      ghost.at(i) = cluster.at(i); // taking initial geometry
      if ( i < poly.size() ){
           //polyhedra too
           putclusterinpoly( cluster.at(i), poly.at(i) );
           Poly ghostpoly = poly.at(i);
           perfectghost( ghostpoly );
           ghost.at(i) = clusterfrompoly( ghostpoly );
      }
	}
}

void Structure::prepimperfectclustersandghosts(){
	
	for ( int i = 0; i < cluster.size(); i++ ){
      formcluster( cluster.at(i) ); // checking
      ghost.at(i) = cluster.at(i); // taking initial geometry
      if ( i < poly.size() ){
           //polyhedra too
           putclusterinpoly( cluster.at(i), poly.at(i) );
           Poly ghostpoly = poly.at(i);
           //perfectghost( ghostpoly );
           //ghost.at(i) = clusterfrompoly( ghostpoly );
      }
	}
}

void Structure:: labelvertices(){
	
  for ( int i = 0 ; i < cluster.size(); i++ ) {
      for ( int j = 0 ; j < cluster.at(i).members.size(); j++ ){
          int k = cluster.at(i).members.at(j);
          //atom k is vertex j in cluster i
          //cerr << "Atom " << k << " is vertex " << j << " in cluster " << i << endl;
          atom.at(k).inclust.push_back( i );
          atom.at(k).isvertex.push_back( j );
      }
  }	

}


void Structure::newcell( vector < double > newcellparam ){
	cell.param = newcellparam;
	cell.makecell();
	cell.makestar();
	
	for ( int i = 0; i < atom.size(); i++){
		atom.at(i).pos = cell.fractocart( atom.at(i).fracpos );
		atom.at(i).initial_pos = atom.at(i).pos;
	}
	
	cell.makegrid();
	fillgrid();
	
	for ( int i = 0 ; i < cluster.size() ; i++ ){
		formcluster( cluster.at(i) ); //updated geometry in cluster
		ghost.at(i).cpos = cluster.at(i).cpos; // update location
		if ( i < poly.size() ){
			putclusterinpoly( cluster.at(i), poly.at(i) ); //update to poly entry
		}
	}
}


//
//newstructure puts a new structure (e.g. xtl waypoint) cell and fracpos into this structure's entities
//it leaves ghosts unchanged
//only call if the new structure is valid, i.e. same number of atoms in clusters
//
void Structure::newstructure( Structure &news ){
	
	cell.param = news.cell.param; // cell overwritten
	cell.makecell();
	cell.makestar();

	for ( int i = 0; i < atom.size(); i++){
		atom.at(i).fracpos = news.atom.at(i).fracpos; // update fractional coords
		atom.at(i).pos = cell.fractocart( atom.at(i).fracpos ); //update cart pos
		atom.at(i).initial_pos = atom.at(i).pos;
	}
	
	cell.makegrid();
	fillgrid();
	
	for ( int i = 0 ; i < cluster.size() ; i++ ){
		formcluster( cluster.at(i) ); //updated geometry in cluster
		ghost.at(i).cpos = cluster.at(i).cpos; // update location
		if ( i < poly.size() ){
			putclusterinpoly( cluster.at(i), poly.at(i) ); //update to poly entry
		}
	}
}

//
//write the GASP polyhedra file with information on ghost and vertex positions
//
void Structure::writepol( string filename ){
     //write out GASP-style polyhedra file for info
     cerr << "Writing structure to " << filename << endl;
     ofstream polfile( filename.c_str() );
     
     polfile << "Polyhedral and cluster position information for " << filename << endl;
     polfile << "CELL VECTORS ijk:" << endl;
     polfile << cell.celli << endl;
     polfile << cell.cellj << endl;
     polfile << cell.cellk << endl;
     polfile << "POLYHEDRON and CLUSTER LIST" << endl;
     
     for ( int i = 0 ; i < cluster.size() ; i++ ){
         polfile << "Poly: " << fixed << setw(6) << i+1 << "   Atoms:";
         for ( int j=0 ; j< cluster.at(i).members.size() ; j++ ){
             polfile << " " << fixed << setw(6) << 1+cluster.at(i).members.at(j);
         }
         polfile << endl; //this completes the line: "Poly: PolID Atoms: catom v1 v2 v3..."
         for ( int j=0 ; j< cluster.at(i).members.size() ; j++ ){
             polfile << fixed << setw(6) << 1+cluster.at(i).members.at(j) << " ";
             polfile << cluster.at(i).cpos + cluster.at(i).bond.at(j) << "   ";
             polfile << ghost.at(i).cpos + ghost.at(i).bond.at(j) << endl;
         }         
     }
     
     polfile << "END POLYHEDRON and CLUSTER LIST" << endl;
     polfile << endl; 
     
     polfile.close();
     return;
}


//
//bond lengths
//
void Structure::writebondlengths( string filename ){
     cerr << "Writing bondlength information to " << filename << endl;
     ofstream bondfile( filename.c_str() );     
     
     for ( int i = 0 ; i < atom.size() ; i++ ){
         int nn = atom.at(i).bondto.size();
         if ( nn < 1 ) continue; //no bonds
         
         for ( int j = 0 ; j < nn ; j++ ){
                 int a = atom.at(i).bondto.at(j);
                 
                 Vector bond1;
                 Vector working;
                 working = atom.at(a).pos - atom.at(i).pos;
                 bond1 = cell.wrapcartdel( working );
                 
                 double dot = sqrt( bond1.sq() ) ;
                 
                 bondfile << i+1 << " " << atom.at(i).species << "  ";
                 bondfile << a+1 << " " << atom.at(a).species << "    ";

                 bondfile << dot << endl;

                 }
     }     
     
     bondfile.close();
}


//
// write all three-body angles
//
void Structure::writeangles( string filename ){
     cerr << "Writing angle information to " << filename << endl;
     ofstream anglefile( filename.c_str() );
     
     for ( int i = 0 ; i < atom.size() ; i++ ){
         int nn = atom.at(i).bondto.size();
         if ( nn < 2 ) continue; //no angle here
         
         for ( int j = 0 ; j < nn ; j++ ){
             for ( int k = j+1 ; k < nn ; k++ ){
                 int a = atom.at(i).bondto.at(j);
                 int b = atom.at(i).bondto.at(k);
                 
                 Vector bond1, bond2;
                 Vector working;
                 working = atom.at(a).pos - atom.at(i).pos;
                 bond1 = cell.wrapcartdel( working );
                 working = atom.at(b).pos - atom.at(i).pos;
                 bond2 = cell.wrapcartdel( working );
                                  
                 bond1.norm(); //using my normer method
                 bond2.norm();
                 
                 double dot = bond1.dot( bond2 );
                 if ( dot > 1.0 ) dot= 1.0;
                 if ( dot < -1.0 ) dot = -1.0;
                 double theta = radtodeg * acos( dot );
                 
                 anglefile << i+1 << " " << atom.at(i).species << "  ";
                 anglefile << a+1 << " " << atom.at(a).species << "  ";
                 anglefile << b+1 << " " << atom.at(b).species << "    ";

                 anglefile << theta << endl;
             }
         }
     }
     
     anglefile.close();
}

void Structure::fitghosts(){
	for ( int i = 0; i < ghost.size(); i++){
		Vector rotor;
		fitrotor( cluster.at(i), ghost.at(i), rotor ); //rotate ghost to fit cluster
	}
}

double Structure::maxrad(){
	double result = 0.;
	for ( int i = 0; i< element.size(); i++ ){
		if ( element.at(i).radius > result ) result = element.at(i).radius;
	}
	return result;
}



//
//contact returns relative positions and contact distances around an atom
//
vector< Contact > Structure::contact( int at, double within ){
        vector< Contact > result;
        result.clear();
        //double lookwithin = 2.0*maxrad;
        Vector diff, deltaf, delta;
        vector< int > near = nearlist( at, within );
        //parse this list for true contacts;
        double r1 = atom.at(at).radius;
        for ( int i = 0 ; i < near.size(); i++ ){
            int ot = near.at(i);
            if ( bonded ( atom.at(at), atom.at(ot) ) ) continue; // not steric
            double r = r1 + atom.at(ot).radius; // contact distance
            
            if ( do14 ){
                 //test isn14 for at and ot
                 if ( isn14( at, ot ) ){
                           //cerr << "Triggering r14scale for atoms " << at+1 << "," << ot+1 << endl;
                           r *= r14scale;
                 }
            }
            
            diff = atom.at(ot).fracpos - atom.at(at).fracpos; // fractional diff;
            deltaf = wrapfracdel( diff ) ; // cell-wrapped;
            delta = cell.fractocart( deltaf ) ;//Cartesian vector to contact ot
            double l2 = delta.sq();
            if ( l2 > r*r ) continue; // not a contact after all
            Contact newcon;
            newcon.dr = delta;
            newcon.l = r;
            newcon.id1 = at;
            newcon.id2 = ot;
            result.push_back( newcon );
            
          
        }
        return result;
}


//
//fullcontact returns relative positions and contact distances for all atoms
//using fullnear list generated in structure...
//
vector< vector< Contact > > Structure::fullcontact ( double within ){
        vector< vector< Contact > > result;
        vector< vector< int > > fullnear = fullnearlist( within );
        result.clear();
        
        int nat = atom.size();
        result.resize( nat );
        Vector diff, deltaf, delta;
        
        for ( int i = 0 ; i < nat ; i++ ){
            double r1 = atom.at(i).radius;
            for ( int j = 0 ; j < fullnear.at(i).size(); j++ ){
                int ot = fullnear.at(i).at(j);
                if ( bonded( atom.at(i), atom.at(ot) ) ) continue; //no steric
                double r = r1 + atom.at(ot).radius;
                
                if ( do14){
                     //test isn14 for ids i and ot
                     if ( isn14( i, ot ) ){
                               //cerr << "Triggering r14scale for atoms " << i+1 << "," << ot+1 << endl;
                               r *= r14scale;
                     }
                }
                
                diff = atom.at(ot).fracpos - atom.at(i).fracpos; // fracdiff
                deltaf = wrapfracdel( diff ) ; // cell-wrapped
                delta = cell.fractocart( deltaf); // Cartesian vector, i to ot
                double l2 = delta.sq();
                if ( l2 > r*r ) continue; // no contact
                Contact newcon;
                newcon.dr = delta;
                newcon.l = r;
                newcon.id1 = i;
                newcon.id2 = ot; // atom notes
                result.at(i).push_back( newcon);
                newcon.dr *= -1; // reverse the polarity
                newcon.id1 = ot;
                newcon.id2 = i; // reverse IDs
                result.at(ot).push_back( newcon );                            
            }
        }

        return result;     
}


//put the bonds from a Bondline object into the atom bondto of a structure
bool Structure::bondsintoatoms( vector < Bondline > &bl ){
	
	for ( int i = 0 ; i < bl.size() ; i++ ){
        int a1 = bl.at( i ).id1;
        int a2 = bl.at( i ).id2;
        int nbars = bl.at(i).bars;
        string sp1 = bl.at(i).ele1;
        string sp2 = bl.at(i).ele2;
         
         //validations!
         if ( a1 >= atom.size() ) {
              cerr << "FATAL: bond input references atom " << a1+1 << " which isn't there." << endl;
              return false;
         }     
         if ( a2 >= atom.size() ) {
              cerr << "FATAL: bond input references atom " << a2+1 << " which isn't there." << endl;
              return false;
         } 
         if ( lowerCaseVersion( sp1 ) != lowerCaseVersion( atom.at(a1).species ) ){
              cerr << "WARNING: bond input disagrees on species of atom " << a1+1;
              cerr << " ( bond " << sp1 << " , xtl " << atom.at(a1).species << " )" << endl;
         }
         if ( lowerCaseVersion( sp2 ) !=  lowerCaseVersion( atom.at(a2).species ) ){
              cerr << "WARNING: bond input disagrees on species of atom " << a2+1;
              cerr << " ( bond " << sp2 << " , xtl " << atom.at(a2).species << " )" << endl;
         }

         //rigidity
         if ( nbars == 6 ){
              atom.at(a1).hybrid = 3; // generic label
              atom.at(a2).hybrid = 3; // generic label
              //cerr << "Locking bond " << a1+1 << " -- " << a2+1 << endl;
         }

         //bonding
         atom.at(a1).bondto.push_back(a2);
         atom.at(a2).bondto.push_back(a1);                  
    }
	return true;
}

void Structure::formbondlines( vector< Bondline > &bl ){
	bl.clear();
	Bondline tempbl;
         
       for ( int i = 0 ; i < atom.size(); i++ ){
           int nb = atom.at(i).bondto.size();
           for ( int j = 0 ; j < nb ; j++ ){
               int foo = atom.at(i).bondto.at(j); // other atom
               if ( foo < i ) continue; // let's not double-report
               tempbl.id1 = i;
               tempbl.id2 = foo;
               // rigidity check!
               tempbl.bars = 5; // default; dihedral is free
			   int hy1 = atom.at(i).hybrid;
			   int hy2 = atom.at(foo).hybrid;	
				
               if ( hy1 != 0 && hy2 != 0 ){
				    if ( hy1 < 0 and hy2 < 0){
						tempbl.bars=5;// no change = special case caught
					}
					else{
						tempbl.bars = 6; // rigid; dihedral locked
					}
               }
               
               tempbl.ele1 = atom.at(i).species;
               tempbl.ele2 = atom.at(foo).species;               
               bl.push_back ( tempbl );
           }           
       }
	
}




void Structure::polybondfind( vector< Polyspec > polyspec, double pad ){
       for ( int i = 0; i < polyspec.size(); i++ ) {
           string csp = lowerCaseVersion( polyspec.at(i).c_species );
           string vsp = lowerCaseVersion( polyspec.at(i).v_species );
           double within = polyspec.at(i).bondlength * pad;
           
           for ( int j = 0 ; j < atom.size(); j++ ){
               if ( lowerCaseVersion( atom.at(j).species ) != csp ) continue; //skip not centre
               atom.at(j).role = -1; // labels poly centre for MC purposes
               vector< int > mine = nearlist( j, within );
               for ( int k = 0; k < mine.size(); k++ ){
                   int foo = mine.at(k);
                   if ( lowerCaseVersion( atom.at(foo).species ) != vsp ) continue;
                   //foo is a potential vertex bonded atom!
                   atom.at(j).bondto.push_back( foo );
                   atom.at( foo ).bondto.push_back( j );
                   //cerr << "Poly-bonding " << j << " to " << foo << endl;
               }
               
           }
           
       }
	
}


void Structure::bondbondfind( vector< Bondspec > bondspec ){
       for ( int i = 0; i < bondspec.size(); i++ ) {
           string asp = lowerCaseVersion( bondspec.at(i).speciesA );
           double within = bondspec.at(i).within;
           
           for ( int j = 0 ; j < atom.size(); j++ ){
               if ( lowerCaseVersion( atom.at(j).species ) != asp ) continue; //skip not centre
               vector< int > mine = nearlist( j, within );
               for ( int k = 0; k < mine.size(); k++ ){
                   int foo = mine.at(k);
                   bool mightbe = false;
                   for ( int baz = 0 ; baz < bondspec.at(i).speciesB.size(); baz++ ){
                       if ( lowerCaseVersion( atom.at(foo).species ) == lowerCaseVersion( bondspec.at(i).speciesB.at(baz) ) ){
                            mightbe = true; //on the list
                       }
                   }
                   if ( !mightbe ) continue;
                   //foo is a potential vertex bonded atom!
                   if ( !isin ( foo, atom.at(j).bondto ) ){
                      atom.at(j).bondto.push_back( foo );
                      atom.at( foo ).bondto.push_back( j );
                      //cerr << "Bonding " << j << " to " << foo << endl;                        
                   }

               }
               
           }
           
       }
}


void Structure::parsepoly(vector< Polyspec > &polyspec ){
	
	
	   for ( int i = 0; i < polyspec.size(); i++ ) {
           polyspec.at(i).nmembers =0 ;//clear to start
           string csp = lowerCaseVersion( polyspec.at(i).c_species );
           string vsp = lowerCaseVersion( polyspec.at(i).v_species );
           
           for ( int j = 0 ; j < atom.size(); j++ ){
               if ( lowerCaseVersion( atom.at(j).species ) != csp ) continue; //skip not centre
               vector< int > mine = atom.at(j).bondto ; // my bonds
               vector< int > vert; // possible vertex atoms
               vert.clear();
               for ( int k = 0; k < mine.size(); k++ ){
                   int foo = mine.at(k);
                   if ( lowerCaseVersion( atom.at(foo).species ) != vsp ) continue; // skip not vertex
                   //cerr << "Atom " << foo << " vertex." << endl;
                   vert.push_back(foo);
               }
               //Have we got a poly?? Verify
               //cerr << "Vertices: " << vert.size() << endl;
               if ( polyspec.at(i).shape == 1 ){
                    //tetrahedron
                    if ( vert.size() != 4 ){
                         cerr << "WARNING: atom " << j+1 << " has " << vert.size() << " vertices, should be 4." << endl;
                    }
                    else{
                         Poly newpoly;
                         newpoly.c_species = csp;
                         newpoly.v_species = vsp;
                         newpoly.shape = 1; //tetrahedron
                         newpoly.bondlength = polyspec.at(i).bondlength;
                         newpoly.members.clear();
                         newpoly.members.push_back( j ); //centre
                         newpoly.cpos = atom.at( j ).pos;
                         newpoly.bond.push_back(nullvec); //will build bonds later!
                         for ( int k = 0 ; k < 4 ; k++ ){
                             newpoly.members.push_back( vert.at( k ) ); // vertices
                             newpoly.bond.push_back(nullvec); // build later
                         }
                         newpoly.type = i ;// this kind
                         polyspec.at( i ).nmembers++ ; // counter
                         poly.push_back ( newpoly );
                         Cluster newcluster = clusterfrompoly ( newpoly );
                         cluster.push_back( newcluster) ; // fill in later
                         ghost.push_back( newcluster) ; // fill in later
                    }
               } // done tet
               if ( polyspec.at(i).shape == 4 ){
                    //square
                    if ( vert.size() != 4 ){
                         cerr << "WARNING: atom " << j+1 << " has " << vert.size() << " vertices, should be 4." << endl;
                    }
                    else{
                         Poly newpoly;
                         newpoly.c_species = csp;
                         newpoly.v_species = vsp;
                         newpoly.shape = 4; //square
                         newpoly.bondlength = polyspec.at(i).bondlength;
                         newpoly.members.clear();
                         newpoly.members.push_back( j ); //centre
                         newpoly.cpos = atom.at( j ).pos;
                         newpoly.bond.push_back(nullvec); //will build bonds later!
                         for ( int k = 0 ; k < 4 ; k++ ){
                             newpoly.members.push_back( vert.at( k ) ); // vertices
                             newpoly.bond.push_back(nullvec); // build later
                         }
                         newpoly.type = i ;// this kind
                         polyspec.at( i ).nmembers++ ; // counter
                         poly.push_back ( newpoly );
                         Cluster newcluster = clusterfrompoly ( newpoly );
                         cluster.push_back( newcluster) ; // fill in later
                         ghost.push_back( newcluster) ; // fill in later
                    }
               } // done square

               if ( polyspec.at(i).shape == 2 ){
                    //octahedron
                    if ( vert.size() != 6 ){
                         cerr << "WARNING: atom " << j+1 << " has " << vert.size() << " vertices, should be 6." << endl;
                    }
                    else{
                         Poly newpoly;
                         newpoly.c_species = csp;
                         newpoly.v_species = vsp;
                         newpoly.shape = 2; //octahedron
                         newpoly.bondlength = polyspec.at(i).bondlength;
                         newpoly.members.clear();
                         newpoly.members.push_back( j ); //centre
                         newpoly.cpos = atom.at( j ).pos;
                         newpoly.bond.push_back(nullvec); //will build bonds later!
                         for ( int k = 0 ; k < 6 ; k++ ){
                             newpoly.members.push_back( vert.at( k ) ); // vertices
                             newpoly.bond.push_back(nullvec); // build later
                         }
                         newpoly.type = i ;// this kind
                         polyspec.at( i ).nmembers++ ; // counter
                         poly.push_back ( newpoly );
                         Cluster newcluster = clusterfrompoly ( newpoly );
                         cluster.push_back( newcluster) ; // fill in later
                         ghost.push_back( newcluster) ; // fill in later
                    }
               } // done oct
               if ( polyspec.at(i).shape == 3 ){
                    //triangle
                    if ( vert.size() != 3 ){
                         cerr << "WARNING: atom " << j+1 << " has " << vert.size() << " vertices, should be 3." << endl;
                    }
                    else{
                         Poly newpoly;
                         newpoly.c_species = csp;
                         newpoly.v_species = vsp;
                         newpoly.shape = 3; //triangle
                         newpoly.bondlength = polyspec.at(i).bondlength;
                         newpoly.members.clear();
                         newpoly.members.push_back( j ); //centre
                         newpoly.cpos = atom.at( j ).pos;
                         newpoly.bond.push_back(nullvec); //will build bonds later!
                         for ( int k = 0 ; k < 3 ; k++ ){
                             newpoly.members.push_back( vert.at( k ) ); // vertices
                             newpoly.bond.push_back(nullvec); // build later
                         }
                         newpoly.type = i ;// this kind
                         polyspec.at( i ).nmembers++ ; // counter
                         poly.push_back ( newpoly );
                         Cluster newcluster = clusterfrompoly ( newpoly );
                         cluster.push_back( newcluster) ; // fill in later
                         ghost.push_back( newcluster) ; // fill in later
                    }
               } // done tri
               if ( polyspec.at(i).shape == 5 ){
                    //bar
                    //special case! make a bar poly out of EACH vert!
                    //cerr << "Working on bars, around atom " << j << endl;
                    //cerr << "Got " << vert.size() << " vertices " << endl;
                    for ( int p =0 ; p < vert.size() ; p++ ){
                         Poly newpoly;
                         newpoly.c_species = csp;
                         newpoly.v_species = vsp;
                         newpoly.shape = 5; //bar
                         newpoly.bondlength = polyspec.at(i).bondlength;
                         newpoly.members.clear();
                         newpoly.members.push_back( j ); //centre
                         newpoly.cpos = atom.at( j ).pos;
                         newpoly.bond.push_back(nullvec); //will build bonds later!
                             newpoly.members.push_back( vert.at( p ) ); // vertices
                             newpoly.bond.push_back(nullvec); // build later
                         newpoly.type = i ;// this kind
                         polyspec.at( i ).nmembers++ ; // counter
                         poly.push_back( newpoly );
                         Cluster newcluster = clusterfrompoly ( newpoly );
                         cluster.push_back( newcluster) ; // fill in later
                         ghost.push_back( newcluster) ; // fill in later
                    }
               } // done bar
               
           }
           
       }
}

vector< Cluster > Structure::candidatecluster ( vector< Bondspec > &bondspec ){

       vector< Cluster > candidate; // initial bonding clusters
       candidate.clear(); 
       
       bool debugme = false;
       
       for ( int i = 0; i < bondspec.size(); i++ ) {
           string asp = lowerCaseVersion( bondspec.at(i).speciesA );
           
           if ( debugme ) cerr << "Checking to bond from species " << asp << endl;
           
           for ( int j = 0 ; j < atom.size(); j++ ){
               if ( lowerCaseVersion( atom.at(j).species ) != asp ) continue; //skip not centre
			   if ( debugme ) cerr << "Checking bonding around atom index "  << j+1 << endl;               
               vector< int > mine = atom.at(j).bondto ; // my bonds
               if ( debugme ){
				   cerr << "Bondto list: ";
				   for ( int q = 0 ; q < mine.size() ; q++ ){
					   cerr << mine.at(q)+1 << " ";
				   }
				   cerr << endl;
			   }
               vector< int > vert; // possible vertex atoms
               vert.clear();
               for ( int k = 0; k < mine.size(); k++ ){
                   int foo = mine.at(k);
                   bool mightbe = false;
                   for ( int baz = 0 ; baz < bondspec.at(i).speciesB.size(); baz++ ){
                       if ( lowerCaseVersion( atom.at(foo).species ) == lowerCaseVersion( bondspec.at(i).speciesB.at(baz) ) ){
                            mightbe = true; //on the list
                            if ( debugme ) cerr << "Good match to atom " << foo+1 << " with species " << atom.at(foo).species << endl;
                       }
                   }

                   if ( !mightbe ) continue;
                   vert.push_back(foo);
               }
               if ( vert.size() == 0 ){
				   if ( debugme ) cerr << "Found no partners." << endl;
				   continue; // I will not make a cluster if I'm all alone!
			   }
               //now we have a cluster - file it
               Cluster newcluster;
               newcluster.members.clear();
               newcluster.members.push_back( j ); // anchor
               newcluster.cpos = atom.at(j).pos;
               newcluster.bond.push_back( nullvec ); // fill in later
               for ( int k = 0 ; k < vert.size() ; k++ ){
                   newcluster.members.push_back( vert.at(k) );
                   newcluster.bond.push_back( nullvec );
               }
               candidate.push_back( newcluster );
               if ( debugme ){
				   cerr << "Built a candidate cluster with members: ";
				   for ( int q = 0 ; q < newcluster.members.size() ; q++ ){
					   cerr << newcluster.members.at(q)+1 << " ";
				   }
				   cerr << endl;
			   }

           }
	
		}
		
	return candidate; // give the list of clusters back to the world.
}

//
//Garibaldi returns a new candidate list of clusters
//by unifying the initial list through adjacent sp2 centres using hybrid label
//
vector< Cluster > Structure::Garibaldi( vector< Cluster > &candidate ){
	vector < Cluster > result; // to return
	
	bool debugme = false;
	
	int nc = candidate.size();
	vector< int > bondto; // unify to me...
	bool gotmates = false; //risorgimento needed?
	
	bondto.resize( nc );
	
	for ( int i = 0; i < nc ; i++ ){
		bondto.at(i) = i; // the identity case
		if ( debugme ) cerr << "Cluster " << i << " gets label " << i << endl;
	}
	
	//start while loop getting bondto
	bool changing = true;
	int rounds = 0;
	
	while ( changing ){
		changing = false;
		rounds++;
		cerr << "Garibaldi performing unification round: " << rounds << endl;
		for ( int i = 0 ; i < nc ; i++ ){
			int c1 = candidate.at(i).members.at(0); // one centre
			if ( atom.at(c1).hybrid == 0 ) continue; //won't unify, not hybrid
			if ( debugme ) cerr << "i-Atom " << c1+1 << " hybrid " << atom.at(c1).hybrid << endl;
			for ( int j = i+1; j < nc ; j++ ){
				int c2 = candidate.at(j).members.at(0);		
				if ( atom.at(c2).hybrid == 0 ) continue; //won't unify, not hybrid	
				if ( !isin( c2 , atom.at(c1).bondto ) ) continue; // I forgot this once. D'oh.
				if ( debugme ) cerr << "j-Atom " << c2+1 << " hybrid " << atom.at(c2).hybrid << endl;			
				if ( atom.at(c1).hybrid == -1 && atom.at(c2).hybrid == -1 ) continue; // special case
				
				//if we're here then some unification is happening...
				gotmates = true;
				//both clusters must be bonded to the lowest available cluster label! which could be either
				int italy = bondto.at(i);
				if ( bondto.at(j) < italy ) italy = bondto.at(j);
				
				if ( italy != bondto.at(i) || italy != bondto.at(j) ) changing = true;
				
				bondto.at(i) = italy;
				bondto.at(j) = italy;
				//so "bondto" labels propagate like a burning algorithm
				if ( debugme ) cerr << "Bonding cluster " << j << " to " << bondto.at(i) << " through " << i << endl;
			}
			if ( debugme ) cerr << "Current bondto for cluster " << i << " is " << bondto.at(i) << endl;
		}
	
	}
	//end while loop getting bondto
	
	
	if ( !gotmates ){
		cerr << "Garibaldi found nothing to do: washing shirts." << endl;
		result = candidate; return result; // don't do anything
	}
	
	//now the unification process must proceed... differently from original GASP
	//do a double pass
	// first pass;all clusters not bonded to self unify with their bondto partner.
	//second pass; only clusters bonded to self are pushed into result
	//we alter candidate in this pass so be careful - shouldn't use it again
	
	for ( int i = 0; i < nc; i++ ){
		if ( bondto.at(i) != i ){
			if ( debugme ) cerr << " Cluster " << i << " unites with " << bondto.at(i) << endl;
			int borg = bondto.at(i);
			Cluster tempclust;
			tempclust = uniclust ( candidate.at(borg), candidate.at(i) );
			candidate.at(borg) = tempclust; // so cluster i has been put into cluster borg
		}
	}

	for ( int i = 0; i < nc; i++ ){
		if ( bondto.at(i) != i ){
			if ( debugme )  cerr << "Not retaining cluster " << i << endl;
			continue;
		} // don't push these;
		if ( debugme ) cerr << "Retaining cluster " << i << endl;
		result.push_back (  candidate.at(i) ); // borg clusters pushed into result
	}
	cerr << "Garibaldi completed risorgimento. Grazie, Giuseppe." << endl;
	return result; // and we're done
}

//
//writemismatches exports the mismatch and clash information
//routine runs a round of poly fitting, mismatch calc and clash finding
//then reports text file
//
void Structure::writemismatches(string filename ){
     cerr << "Calculating mismatch and clash info;" << endl;
     cerr << "writing to file " << filename << endl;
     
     ofstream misout( filename.c_str() );
     stringstream misstream; // header summary info
     misstream << "Measure: |mismatch| mismatch^2 |clash| clash^2" << endl;
     
     stringstream clashstream; // clash list
     
     //fit the ghost positions to the current clusters
	 fitghosts(); // all done by structure method
	 cerr << "Ghosts fitted." << endl;
     //obtain mismatches, grading worst magnitude and square, and totalling
     double worstmodm =0.0; //mismatch length
     double worstm2 = 0.0; //squared mismatch (penalty)
     double worstc = 0.0; //clash overlap distance
     double worstc2 = 0.0; //overlap squared (penalty)
     double totalmodm = 0.0;
     double totalm2 = 0.0;
     double totalc = 0.0;
     double totalc2 = 0.0;
     int nmis = 0;
     int nc = 0;
     
     for ( int i = 0 ; i < atom.size(); i++ ){
         Vector mis;
         int ncl = atom.at(i).inclust.size();
         for ( int wc = 0 ; wc < ncl ; wc++ ){
             int ac = atom.at(i).inclust.at(wc);
             int av = atom.at(i).isvertex.at(wc);
             mis = mymis ( atom.at(i), ghost.at(ac), av );
             nmis++; // running count of how many mismatches exist in structure!
             //got mismatch vector. Check size;
             double mis2 = mis.sq();
             double modmis = sqrt( mis2 );
             totalmodm += modmis; // running count
             totalm2 += mis2;
             
             if ( mis2 > worstm2 ){
                  worstm2 = mis2;
                  worstmodm = modmis;
             }
         }
     }
     cerr << "Mismatches measured" << endl;
     
     //get clash information
     double lookie = 2.0*maxrad(); // checks elements for max radius
     vector< vector< Contact > > thecon = fullcontact( lookie ); // calls fullnearlist 
     cerr << "Clashes identified" << endl; 
     
     for ( int i = 0 ; i < atom.size(); i++ ){
         vector< Contact > cont = thecon.at(i); // this atom's contacts
         int ncon = cont.size();
         for ( int wcon = 0; wcon < ncon ; wcon++ ){
             int ot = cont.at(wcon).id2; // the other atom
             if ( ot < i ) continue; //no need to doublecount
             nc++;
             
             double r = sqrt( cont.at(wcon).dr.sq() ); // true distance
             double c = cont.at(wcon).l - r ;// clash overlap distance
             clashstream << i+1 << " " << atom.at(i).species << " ";
             clashstream << ot+1 << " " << atom.at(ot).species << " ";
             clashstream << c << " " << r << " " << cont.at(wcon).l << endl;
             //written a clash line to clashstream
             
             totalc += c; //running count
             totalc2 += c*c;
             if ( c > worstc ){
                  worstc = c;
                  worstc2 = c*c;
             }        
             
         }
         
     }     
     //write header stream
     misstream << "Worst " << worstmodm << " " << worstm2 << " ";
     misstream << worstc << " " << worstc2 << endl;
     
     misstream << "Total " << totalmodm << " " << totalm2 << " ";
     misstream << totalc << " " << totalc2 << endl;
     
     misstream << "N " << nmis << " " << nmis << " ";
     misstream << nc << " " << nc << endl;
     
     //header to file
     misout << misstream.str();
     //clashes to file
     misout << "Clash list:" << endl;
     misout << clashstream.str();
                         
    misout.close();
    return;
}

void Structure::jiggle( double stepsize ){
	cerr << "Performing random jiggle before relaxation." << endl;
	for ( int i = 0; i < atom.size(); i++ ){
		Vector jigglemove = unit_vector()*stepsize; //boing
		
		Vector newpos = atom.at(i).pos + jigglemove;
		atom.at(i).pos = cell.wrapcartpos( newpos );
		atom.at(i).fracpos = cell.carttofrac( atom.at(i).pos );
	}
	for ( int i = 0; i < cluster.size(); i++ ){
		formcluster( cluster.at(i) );
		ghost.at(i).cpos = cluster.at(i).cpos;
	}
	cerr << "Jiggling process complete." << endl;
}


//
//update does atom and cluster positions based on move vVector
//
void Structure::update( vector< Vector > &move ){
	for ( int i = 0; i < atom.size(); i++){
		Vector newpos = atom.at(i).pos + move.at(i); // raw new position
		atom.at(i).pos = cell.wrapcartpos( newpos ); // cell pbc handling
		atom.at(i).fracpos = cell.carttofrac( atom.at(i).pos );
	}
	for ( int i = 0; i < cluster.size(); i++){
		formcluster( cluster.at(i) ); // update cluster geometry
		ghost.at(i).cpos = cluster.at(i).cpos; // ghost location
	}
	
}   


void Structure::shift( Vector move ){
	for ( int i = 0; i < atom.size(); i++){
		cerr << i << " " << atom.at(i).fracpos << endl;
		Vector newpos = atom.at(i).fracpos + move ; 
		atom.at(i).fracpos = wrapfracpos( newpos ); // cell pbc handling
		cerr << i << " " << atom.at(i).fracpos << endl << endl;
	}	
}


bool Structure::read_coords( string filename){

  cerr << "Reading structure from " << filename << endl;

  ifstream coo_file(filename.c_str() );

  string line_in;
  vector< string > read_in;

  atom.clear();
  cell.param.clear();
  cell.param.resize(6);

  Vector fracpos;
  Atom dummyatom;
  int nlines=0;

  //flags to mark stages:
  bool startread=true;
  bool abcread = false;
  bool angleread = false;
  bool midread = false;
  bool atomread = false;
  int dashlines = 0;

  while( getline(coo_file, line_in) ){
    nlines++;
    dummyatom = nullatom;
    cerr << "Reading line " << nlines << ": ";
    cerr << line_in << endl; // dump to screen as we go
    read_in = chop(line_in, " ");

    if ( atomread ){
         if ( read_in.at(0).substr(0,1) == "-" ){
              //dash line
              dashlines++;
              cerr << "Found dashed line " << dashlines << endl;
              if ( dashlines == 2 ){
                   cerr << "That is the end of the atom information block " << endl;
                   break; // and leave the read loop
              }
              continue; // go on to next line...
         }
         cerr << "Atom information line: " << endl;
         dummyatom.species = to_lower(read_in[0]);
         // skip the label column
         fracpos.x = atof(read_in[2].c_str());  
         fracpos.y = atof(read_in[3].c_str());  
         fracpos.z = atof(read_in[4].c_str());  

	 Vector inpos = wrapfracpos( fracpos );
         dummyatom.fracpos = inpos;
         atom.push_back(dummyatom);         
    }
    else if ( midread ){
         if ( read_in.at(0) == "Elmt" ){
            cerr << "Reached atom information block" << endl;
            atomread = true; //ready to read next line     
         }
    }
    else if ( angleread ){
         if ( read_in.size() != 10 ){
            cerr << "Misshapen angles line? I die now." << endl;
            return( false );
         }
         string stringalpha = read_in.at(2);
         string stringbeta = read_in.at(5);
         string stringgamma = read_in.at(8); // warning these have semicolons on them
         
         stringalpha.erase( stringalpha.end()-1, stringalpha.end() );
         stringbeta.erase( stringbeta.end()-1, stringbeta.end() );
         stringgamma.erase( stringgamma.end()-1, stringgamma.end() );
         cerr << "Trimmed angles: " << stringalpha << " " << stringbeta << " " << stringgamma << endl;
         cell.param.at(3) = atof ( stringalpha.c_str() );
         cell.param.at(4) = atof ( stringbeta.c_str() );
         cell.param.at(5) = atof ( stringgamma.c_str() );         
         
         midread = true; //ready to skip to atom information                
    }
    else if ( abcread ){
         if ( read_in.size() != 10 ){
            cerr << "Misshapen a b c line? I die now." << endl;
            return false;;
         }
         string stringa = read_in.at(2);
         string stringb = read_in.at(5);
         string stringc = read_in.at(8); // warning these have semicolons on them
         
         stringa.erase( stringa.end()-1, stringa.end() );
         stringb.erase( stringb.end()-1, stringb.end() );
         stringc.erase( stringc.end()-1, stringc.end() );
         cerr << "Trimmed a b c: " << stringa << " " << stringb << " " << stringc << endl;
         cell.param.at(0) = atof ( stringa.c_str() );
         cell.param.at(1) = atof ( stringb.c_str() );
         cell.param.at(2) = atof ( stringc.c_str() );
         
         angleread = true; //ready to read next line    
    }
    else if ( startread ){
         if ( read_in.at(0) == "Unit" ){
            cerr << "Reached Unit cell parameters" << endl;
            abcread = true; //ready to read next line     
         }
    }
    else{
         cerr << "How the hell did you get here?" << endl;
         return( false );
    }


  }

  cerr << "Done reading structure." << endl;
  cerr << "Found " << atom.size() << " atoms." << endl;

  coo_file.close();
  
  return true;
} // ENDOFREAD_coords


void Structure::stripdupes( double prec ){
        //prec defaults to 0.0005 by definition in header
	vector< bool > isduped;
	isduped.resize( atom.size() ); // set of labels
	
	vector< Atom > newatom;
	newatom.clear();
	
	for ( int i = 0; i < atom.size(); i++ ){
		isduped.at(i) = false; // initialise
	}
	for ( int i = 0; i < atom.size(); i++ ){
		if ( isduped.at(i) ) continue; // already marked for death
		for ( int j = i+1; j < atom.size(); j++ ){
			//double eps = 5.E-4; // tiny fractional difference
                        double eps = prec; //tolerance on frac difference
			Vector delta; // coordinate difference
			delta = atom.at(i).fracpos - atom.at(j).fracpos;
			Vector wdel = wrapfracdel( delta ); // trap for across cell boundary
			double mydot = wdel.sq();
			if ( mydot < eps ){
				cerr << "Atom " << j+1 << " has same position as atom " << i+1 << endl;
				isduped.at(j) = true;
			}
		}
	}
	//copy non-duped atoms into newatom vector
	for ( int i = 0 ; i <  atom.size()  ; i++ ){ 
		if ( isduped.at(i) ) continue; // already checked
		newatom.push_back( atom.at(i) );
	}	
	atom = newatom; // replace old coord list with new.
}

void Structure::super( int nx, int ny, int nz ){

  vector< Atom > newatom;
  newatom.clear();
  for ( int i = 0 ; i < atom.size() ; i++ ){
      for ( int mx = 0 ; mx < nx; mx++ ){
          for ( int my = 0 ; my < ny ; my++ ){
              for ( int mz = 0 ; mz < nz ; mz++ ){
                  Atom copyatom = atom.at(i);
                  copyatom.fracpos.x += 1.0*mx; //copy up
                  copyatom.fracpos.x /= nx; //rescale
                  copyatom.fracpos.y += 1.0*my; //copy up
                  copyatom.fracpos.y /= ny; //rescale
                  copyatom.fracpos.z += 1.0*mz; //copy up
                  copyatom.fracpos.z /= nz; //rescale
                  newatom.push_back(copyatom ) ; // on the list
              }
          }
      }
  }
  
  atom = newatom; // replace with the new list
  cell.param.at(0) *= nx;
  cell.param.at(1) *= ny;
  cell.param.at(2) *= nz; // remember to actually expand cell!
  	
}
void Structure::othersuper( int nx, int ny, int nz ){

  vector< Atom > newatom;
  newatom.clear();

      for ( int mx = 0 ; mx < nx; mx++ ){
          for ( int my = 0 ; my < ny ; my++ ){
              for ( int mz = 0 ; mz < nz ; mz++ ){
                for ( int i = 0 ; i < atom.size() ; i++ ){                 
                  Atom copyatom = atom.at(i);
                  copyatom.fracpos.x += 1.0*mx; //copy up
                  copyatom.fracpos.x /= nx; //rescale
                  copyatom.fracpos.y += 1.0*my; //copy up
                  copyatom.fracpos.y /= ny; //rescale
                  copyatom.fracpos.z += 1.0*mz; //copy up
                  copyatom.fracpos.z /= nz; //rescale
                  newatom.push_back(copyatom ) ; // on the list
              }
          }
      }
  }
  
  atom = newatom; // replace with the new list
  cell.param.at(0) *= nx;
  cell.param.at(1) *= ny;
  cell.param.at(2) *= nz; // remember to actually expand cell!
  	
}
   

void Structure::writecif(string filename, string title){
	cerr << "Writing structure to " << filename << endl;
  ofstream cif_out(filename.c_str());
  cif_out << "# " << title << endl << endl;
  cif_out << "data_" << filename << endl << endl;
  
  cif_out << "_cell_length_a          " << cell.param.at(0) << endl;
  cif_out << "_cell_length_b          " << cell.param.at(1) << endl;  
  cif_out << "_cell_length_c          " << cell.param.at(2) << endl;
  cif_out << "_cell_angle_alpha       " << cell.param.at(3) << endl;
  cif_out << "_cell_angle_beta        " << cell.param.at(4) << endl;
  cif_out << "_cell_angle_gamma       " << cell.param.at(5) << endl;
  cif_out << "_symmetry_space_group_name_H-M 'P 1'" << endl;
  cif_out << "_symmetry_Int_Tables_number 1" << endl << endl;

  cif_out << "loop_" << endl;
  cif_out << "_atom_site_type_symbol" << endl;
  cif_out << "_atom_site_label" << endl;
  cif_out << "_atom_site_fract_x" << endl;
  cif_out << "_atom_site_fract_y" << endl;
  cif_out << "_atom_site_fract_z" << endl;

  //loop now
  int sitelabel = 0;
  for ( int acount = 0 ; acount < atom.size(); acount++ ){
      sitelabel = acount+1;
      cif_out << atom.at(acount).species << " ";
      cif_out << atom.at(acount).species << "_" << sitelabel << " ";
      cif_out << atom.at(acount).fracpos << endl;
  }

  cif_out.close();
} //ENDOFWRITEcif

bool Structure::readcif(string filename ){
	 cerr << "Seeking to read structure from " << filename << endl;
     atom.clear();
     cell.param.clear();
     cell.param.resize(6);
     
     Vector fracPos;
     Atom dummyAtom;
     //bool reading_atoms = false;
     bool hasCell = false;
     bool hasAtoms = false;
     
     string lineIn;
     string foo; //use as buffer
     vector< string > tokens;
     vector< string > linesIn;
     vector< vector< string > > tokensArray;
     vector< int > lineStatus; //label content by type
     
	//new logic:
	//first, read lines into an array of lines without processing
    ifstream inFile( filename.c_str() );
    if ( inFile.fail() ) {
       cerr << "Problem! Cannot read file " << filename << endl;
       return false;     
    }

    while ( getline(inFile, lineIn)  ) {
		linesIn.push_back(lineIn); //filling the array of lines
		lineStatus.push_back(0); //blank status token
	}
	inFile.close();
	
	//second, clean and tokenize non-blank lines into array of arrays of tokens
	for ( int i = 0 ; i < linesIn.size() ; i++ ){
		//cerr << linesIn.at(i) << " ::: " << linesIn.at(i).length() << endl;
		if ( linesIn.at(i).length() == 0 ){
			lineStatus.at(i) = 99; //true blank
			vector<string> dummy;
			tokensArray.push_back(dummy);
		}
		else{
			tokens = chop( linesIn.at(i), " \t" );
			tokensArray.push_back( tokens );
		}
	}
	
	//third, identify cell parameter lines, loop starts, atoms site info lines
	for ( int i = 0 ; i < tokensArray.size() ; i++){
		if ( tokensArray.at(i).size() == 0 ){
			lineStatus.at(i) = 99; //whitespace line
			continue; //is blank and already labelled
		}
		string entry1 = lowerCaseVersion( tokensArray.at(i).at(0) );
		if (  entry1.find("loop_") != string::npos ){
			lineStatus.at(i) = 1; //loop start line
			cerr << "Loop start at line " << i << endl;
		}
		else if (  entry1.find("_cell_length") != string::npos ){
			lineStatus.at(i) = 2; //cell length
			cerr << "Cell length at line " << i << endl;
			if (tokensArray.at(i).size() < 2){
				cerr << "ERROR: cell length label without value." << endl;
				return false;
			}
			else{
				if ( entry1 == "_cell_length_a" ){
					double ca = atof( tokensArray.at(i).at(1).c_str() );
					cell.param.at(0) = ca;
				}
				else if ( entry1 == "_cell_length_b" ){
					double cb = atof( tokensArray.at(i).at(1).c_str() );
					cell.param.at(1) = cb;
				}
				else if ( entry1 == "_cell_length_c" ){
					double cc = atof( tokensArray.at(i).at(1).c_str() );
					cell.param.at(2) = cc;
				}
				else {
					cerr << "ERROR: cell length not a,b, or c. That shouldn't happen." << endl;
					return false;
				}
			}
		} 
		else if (  entry1.find("_cell_angle") != string::npos ){
			lineStatus.at(i) = 3; //cell angle
			cerr << "Cell angle at line " << i << endl;
			if (tokensArray.at(i).size() < 2){
				cerr << "ERROR: cell angle label without value." << endl;
				return false;
			}
			else{
				if ( entry1 == "_cell_angle_alpha" ){
					double aa = atof( tokensArray.at(i).at(1).c_str() );
					cell.param.at(3) = aa;
				}
				else if ( entry1 == "_cell_angle_beta" ){
					double ab = atof( tokensArray.at(i).at(1).c_str() );
					cell.param.at(4) = ab;
				}
				else if ( entry1 == "_cell_angle_gamma" ){
					double ac = atof( tokensArray.at(i).at(1).c_str() );
					cell.param.at(5) = ac;
				}
				else {
					cerr << "ERROR: cell angle not alpha, beta or gamma. That shouldn't happen." << endl;
					return false;
				}
			}
		}
		else if (  entry1.find("_atom_site") != string::npos ){
			lineStatus.at(i) = 4; //atom site specifier
			cerr << "Atom site info at line " << i << endl;
		}
	}
	
	//fourth, parse info
	for ( int i = 0 ; i < tokensArray.size() ; i++ ){
		if ( lineStatus.at(i) != 1 ){
			continue; //not a loop start
		}
		cerr << "Checking loop starting at line " << i << endl;
		int siteCount = 0;
		int xTag = -1;
		int yTag = -1;
		int zTag = -1;
		//int qTag = -1;
		int elTag = -1;
		for ( int j = i+1; j < tokensArray.size() ; j++ ){
			if ( lineStatus.at(j) == 1 || lineStatus.at(j) == 99 ){
				cerr << "Ending read." << endl;
				break; //next loop terminates j loop, as does blank line
			}
			else if ( lineStatus.at(j) == 2 || lineStatus.at(j) == 3 ){
				cerr << "Ending read." << endl;
				break; //cell info lines shouldn't be in an atom loop
			}
			//if we got here we should be either reading:
			//atom site info tag; atom line; something broken
			if ( lineStatus.at(j) == 4 ){
				siteCount++; //increment the site count for this loop
				string voo = lowerCaseVersion( tokensArray.at(j).at(0) );
				if ( voo == "_atom_site_type_symbol" ){
					cerr << "Type symbol recognised as item " << siteCount << " in atom line." << endl;
					elTag = siteCount - 1;
				}
				else if ( voo == "_atom_site_fract_x" ){
					cerr << "X fractional coordinate recognised as item " << siteCount << " in atom line." << endl;
					xTag = siteCount - 1;							
				}
				else if ( voo == "_atom_site_fract_y" ){
					cerr << "Y fractional coordinate recognised as item " << siteCount << " in atom line." << endl;
					yTag = siteCount - 1;							
				}
				else if ( voo == "_atom_site_fract_z" ){
					cerr << "Z fractional coordinate recognised as item " << siteCount << " in atom line." << endl;
					zTag = siteCount - 1;							
				}
				//else if ( voo == "_atom_site_charge" ){
					//cerr << "Charge recognised as item " << siteCount << " in atom line." << endl;
					//qTag = siteCount - 1;							
				//}				
			}
			else if ( tokensArray.at(j).size() != siteCount ){
				cerr << "Line " << j << " ::: " << linesIn.at(j) << endl;
				cerr << "Wrong number of entries for an atom line? Ending read." << endl;
				break; //done with atom loop
			}
			else {
				//should be a good atom line!
				dummyAtom = nullatom;
				if ( elTag > -1 ) dummyAtom.species = tokensArray.at(j).at( elTag ); // read species
				////if ( qTag > -1 ) dummyAtom.charge = atof( tokensArray.at(j).at( qTag ).c_str() ); //read charge if we have one
				if ( xTag > -1 ) dummyAtom.fracpos.x =  atof( tokensArray.at(j).at( xTag ).c_str() ) ; 
				if ( yTag > -1 ) dummyAtom.fracpos.y =  atof( tokensArray.at(j).at( yTag ).c_str() ) ; 
				if ( zTag > -1 ) dummyAtom.fracpos.z =  atof( tokensArray.at(j).at( zTag ).c_str() ) ;
				
				//cerr << "Atom line properties: sp " << dummyAtom.species;
				//cerr << " xyz " << dummyAtom.fracpos;
				//cerr << " q " << dummyAtom.charge << endl;
				Vector conv = dummyAtom.fracpos;					
				dummyAtom.fracpos = wrapfracpos( conv );
				atom.push_back( dummyAtom ); //add to array!				
			}
			
		}
	}

    cerr << "CELL: " << cell.printcell() << endl;
    hasCell = true;
    for ( int j = 0; j < 6; j++ ){
		if ( almostEqual( cell.param.at(j), 0. ) ){
			hasCell=false;
		}
	}	
	if ( atom.size() > 0 ) hasAtoms = true;
	cerr << "Found " << atom.size() << " atom entries." << endl; 

	if ( hasCell && hasAtoms ){
		return true; // all good
	}
	return false;	
	
}


//
//write the MCGasp cluster file with elements and vertex positions
//
void Structure::writeclusters( string filename ){
     //write out GASP-style polyhedra file for info
     cerr << "Writing molecular clusters to " << filename << endl;
     ofstream polfile( filename.c_str() );
     
     polfile << "Molecular cluster information for " << filename << endl;
     
     for ( int i = 0 ; i < cluster.size() ; i++ ){
		 polfile << "CLUSTER with atoms: ";
         for ( int j=0 ; j< cluster.at(i).members.size() ; j++ ){
             polfile << " " << fixed << setw(6) << 1+cluster.at(i).members.at(j);
         }
         polfile << endl; //this completes the line: "CLUSTER with atoms: id1 id2 id3 ..."
         for ( int j=0 ; j< cluster.at(i).members.size() ; j++ ){
			 //write the element
			 string ele =atom.at(  cluster.at(i).members.at(j) ).species;
			 polfile << ele << " ";
			 //write the Cartesian position, not wrapped (local cluster)
             polfile << cluster.at(i).cpos + cluster.at(i).bond.at(j) << endl;
         }
     }
     polfile << "END";
     polfile << endl; 
     polfile.close();
     return;
}
  
Vector Structure::randompos(){
	Vector randomf = cubic_vector();
	Vector randomw = wrapfracpos( randomf ); // wrapped into cell;
	Vector result = cell.fractocart( randomw ); // made Cartesian
	return result; // a random Cartesian position in box
}

void Structure::find14(){
     //for each atom
     //neighbours are in bondto
     //these are 1-2
     //build 1-3 arrays:
     //these contain bondto neighbours of 1-2 which are not also 1-2 or self
     //then build 1-4 arrays:
     //these contain bondto neighbours of 1-3 which are not also 1-2 or 1-3 or self
     //these are placed in the n14 arrays of both atoms if they are not already in there.
     //use the isin bool util :)
     vector< int > n12t;
     vector< int > n13t;
     vector< int > n14t;

     for ( int i =0; i < atom.size() ; i++ ){
         atom.at(i).n14.clear();
     }
     
     for ( int i =0; i < atom.size() ; i++ ){
         n12t.clear();
         n13t.clear();
         n14t.clear();
         
         n12t = atom.at(i).bondto; //copy for clarity!
         
         int n2 = n12t.size();
         for ( int j = 0; j < n2 ; j++ ){
             int o2 = n12t.at(j); //a neighbour atom
             int o2n = atom.at(o2).bondto.size();
             
             for ( int k = 0 ; k < o2n ; k++ ){
                 int o3 = atom.at(o2).bondto.at(k); //potentially a 1-3
                 if ( isin ( o3, n12t ) ) continue; //is already n12, pass
                 if ( o3 == i ) continue; //self! pass.
                 if ( isin( o3, n13t ) ) continue; // already listed
                 n13t.push_back( o3);
             }
         }
         
         //n13t is now fully populated.
         //repeat search seeking n14.
         int n3 = n13t.size();
         for ( int j = 0; j < n3 ; j++ ){
             int o3 = n13t.at(j); //a neighbour atom
             int o3n = atom.at(o3).bondto.size();
             
             for ( int k = 0 ; k < o3n ; k++ ){
                 int o4 = atom.at(o3).bondto.at(k); //potentially a 1-3
                 if ( isin( o4, n12t ) ) continue; //is already n12, pass
                 if ( isin ( o4, n13t ) ) continue; //is already n13, pass
                 if ( o4 == i ) continue; //self! pass. Shouldn't happen?
                 if ( isin( o4, n14t ) ) continue; // already listed
                 n14t.push_back( o4);
             }
         }         
         
         //REPORT step for debugging
         cerr << "1-4 NEIGHBOURS of atom " << i+1 << ": ";
         for ( int j = 0; j < n14t.size() ; j++ ){
             cerr << " " << 1+ n14t.at(j);
         }
         cerr << endl;
         
         //now pass info from n14 into atoms
         for ( int j = 0 ; j < n14t.size() ; j++ ){
             int o4 = n14t.at( j ); //other atom to i
             if ( !isin( o4, atom.at(i).n14 ) ){
                  atom.at(i).n14.push_back( o4 );
             }
             if ( !isin( i, atom.at(o4).n14 ) ){
                  atom.at(o4).n14.push_back( i );
             }
         }
         
     }
            
            
     
}


bool Structure::isn14( int id1, int id2 ){
     //return true if id2 is in the n14 array of atom(id2)
     if ( isin( id2, atom.at(id1).n14 ) ) return true;
     return false;
}


