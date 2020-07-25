#include "Rotors.h"

vector< Cluster > clus1; // real position clusters
vector< Cluster > ghost1; //ideal position clusters
Vector rotor;
vector< double > mbefore; // fitting mismatches
vector< double > mstretch, mbend; //bond length change effect, residual (bend) effect

bool readpol( string filename, vector< Cluster > &clus, vector< Cluster > &ghost ); // returns false if the read is bad.

int main(int argc, char *argv[]){
    
  string infilename;
  
  //cerr << "Args: " << argc << endl;
  if ( argc < 3 ) {
     stringstream whatn;
     whatn << argv[0];
     string myname = whatn.str();
     cerr << "Insufficient filenames given." << endl;
     cerr << "Please give a filename for a pol file for readin," << endl;
     cerr << " and a filename for the analysis output." << endl;
     cerr << " thus: " << myname << " infile outfile" << endl;
     cerr << " Will report cor each polyhedron: " << endl;
     cerr << " the mismatch between its real and ghost vertex positions," << endl;
     cerr << " and the bending and stretching components of that mismatch." << endl;
     exit(0); 
  }
  
  stringstream what1, what2;
  what1 << argv[1];
  what2 << argv[2];
  
  string name1, name2;
  name1 = what1.str();
  name2 = what2.str();

  //bool goodread = readpol ( infilename, clus1, ghost1 );
  bool good1 = readpol ( name1, clus1, ghost1 );
  
  if (!good1){
     cerr << "Problem with reading " << name1 << endl;
     exit(1);  
  }    

  int n1 = clus1.size();

  //loop and solve cluster by cluster ...
  for ( int i = 0 ; i < n1 ; i++ ){
      int m1 = clus1.at(i).bond.size();
 
      // first do vector mismatch before rotation fitting;
      double mis = 0.0;
      Vector diff;
      
      for (int j = 0 ; j < m1 ; j++ ){
          diff = clus1.at(i).bond.at(j) - ghost1.at(i).bond.at(j);
          mis += diff.sq();
          
      }      
      cerr << i+1 << " : Mismatch-squared between real and ghost: " << mis << endl;
      mbefore.push_back( mis );
      
      //and finally do the stretch versus bend breakdown.
      double stretches = 0.;
      for (int j = 0 ; j < m1 ; j++ ){
		  double l1 = sqrt( clus1.at(i).bond.at(j).sq() ); //bond length 1
		  double l2 = sqrt( ghost1.at(i).bond.at(j).sq() ); //bond length 2
		  stretches += ( l1-l2)*(l1-l2); //square of change in bond length          
      }       
      
      cerr << i+1 << " : Stretching mismatch (bond-length-change component): " << stretches << endl;
      mstretch.push_back( stretches );
      double residual = mis - stretches;
      cerr << i+1 << " : Bending mismatch (non-bond-length-change component): " << residual << endl;      
      mbend.push_back( residual );
      
      
  }

  //write to output file
  ofstream outfile( name2.c_str() );
  outfile << "ID M2-total M2-bend M2-stretch" << endl;
  for ( int i = 0 ; i < n1 ; i++ ){
      outfile << i+1 << " " << mbefore.at(i) << " ";
      outfile << mbend.at(i) << " ";
      outfile << mstretch.at(i) << endl;
  }  
  
  outfile.close();
        
  exit(0); //final exit
} // END OF MAIN


//
//readpol pulls clusters data out of the given file
//
bool readpol( string filename, vector< Cluster > &clus, vector< Cluster > &ghost ){
     //read initial_structure vector< Atom > and initial_cellparam from input_structure_name
     cerr << "Seeking to read structure from " << filename << endl;

     clus.clear();
     ghost.clear();
     
     Vector fracpos;
     bool reading_clusters = false;
     
     
     string linein;
     string foo; //use as buffer
     vector< string > tokensin;
     
     vector< vector< string > > clustertoken; // dump everything in here!
     
     ifstream polfile( filename.c_str() );
     if ( polfile.fail() ) {
        cerr << "Problem! Cannot read file " << filename << endl;
        return false;     
     }

     while( getline(polfile, foo) ){
            linein = to_lower(foo);// lose case for simplicity
            //cerr << linein << endl;
            tokensin = chop(linein, " \t");
            if (tokensin.size() == 0 ) continue; //blank line?
            
            if (tokensin.at(0) == "end" ) {
               cerr << "END of pol file" << endl;
               break; // reached end of xtl format
            }
            else if ( tokensin.at(0) == "polyhedron" ) {
                 reading_clusters = true;
                 cerr << "START of cluster list" << endl;
                 continue; //
            }            
            else if ( reading_clusters ) {
                 clustertoken.push_back(tokensin); //save the line for processing
                 //cerr << linein << endl;
            }
            
     } // end of while
    
    polfile.close();
    cerr << "Read pol file into array." << endl;
    //okay, got text into array. Process:
    int entry = 0;
    int lines = clustertoken.size();
    bool working = true;

    Cluster dummyclus, dummyghost;

    while (working) {
          entry++;
          //cerr << "Entry " << entry << endl;
          if ( entry > lines ){
               working = false;
               cerr << "End of data." << endl;
               break; // leave reading
          }
          string key = clustertoken.at(entry-1).at(0);
          if ( key == "poly:" ){

                    dummyclus = nullcluster; // blank entry
                    dummyghost = nullcluster;
                    int id = atoi( clustertoken.at(entry-1).at(1).c_str());
                    //cerr << "Reading cluster " << id << endl;
                    int len = clustertoken.at(entry-1).size();
                    len -= 3 ;// number of atom entries
                    //cerr << len << " members." << endl;
                    for ( int j =3 ; j < len+3 ; j++ ){
                        id = atoi( clustertoken.at(entry-1).at(j).c_str() ); //offset not needed
                        //cerr << "Member id " << id << endl;
                        dummyclus.members.push_back( id );
                        dummyghost.members.push_back( id );
                    }
                    //now loop reading lines
                    for ( int j = 0 ; j < len ; j++ ){
                        entry++;
                        //cerr << "Entry " << entry << endl;
                        if ( entry > lines ){
                             cerr << "Fallen off the edge of the world." << endl;
                             return false;
                        }
                        if ( clustertoken.at(entry-1).size() < 7 ){
                             cerr << "Malformed atom line in cluster." << endl;
                             return false;
                        }
                        if ( clustertoken.at(entry-1).at(0) == "poly:" ){
                             cerr << "Problem: hit Poly: line during read of previous entry!" << endl;
                             return false;
                        }
                        id = atoi( clustertoken.at(entry-1).at(0).c_str() ); // should be atom id
                        //cerr << "ID " << id << ": ";
                        if ( id != dummyclus.members.at(j) ){
                             cerr << "Mismatched atoms? ";
                             cerr << id << " v. " << dummyclus.members.at(j) << endl;
                             return false;
                        }
                        double xc, yc, zc, xg, yg, zg;
                        xc = atof( clustertoken.at(entry-1).at(1).c_str() );
                        yc = atof( clustertoken.at(entry-1).at(2).c_str() );
                        zc = atof( clustertoken.at(entry-1).at(3).c_str() ); 
                        xg = atof( clustertoken.at(entry-1).at(4).c_str() );                                              
                        yg = atof( clustertoken.at(entry-1).at(5).c_str() );
                        zg = atof( clustertoken.at(entry-1).at(6).c_str() );
                        Vector cpos = Vector( xc, yc, zc );
                        Vector gpos = Vector( xg, yg, zg );
                        dummyclus.bond.push_back( cpos );
                        dummyghost.bond.push_back( gpos ); 
                        //cerr << cpos << " / " << gpos << endl;               
                    }
          clus.push_back( dummyclus );
          ghost.push_back( dummyghost );          
          }
          
          
          
    } // done working on parse

    cerr << "Found " << clus.size() << " polyhedra." << endl;
    //do cpos generation here
    int np = clus.size();
    for ( int i = 0 ; i < np ; i++ ){
        Vector newc = nullvec;
        int nm = clus.at(i).bond.size();
        for ( int j = 0 ; j < nm ; j++ ){
            newc += clus.at(i).bond.at(j);
        }
        newc /= nm; // got average position
        for ( int j = 0 ; j < nm ; j++ ){
            clus.at(i).bond.at(j) -= newc; // make relative bond
        }
        clus.at(i).cpos = newc;
        newc = nullvec; // re-zero
        for ( int j = 0 ; j < nm ; j++ ){
            newc += ghost.at(i).bond.at(j);
        }
        newc /= nm; // got average position
        for ( int j = 0 ; j < nm ; j++ ){
            ghost.at(i).bond.at(j) -= newc; // make relative bond
        }
        ghost.at(i).cpos = newc;
    }
    

    return true;
}
