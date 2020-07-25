#include "Rotors.h"

vector< Cluster > clus1, clus2; // real position clusters
vector< Cluster > ghost1, ghost2; //ideal position clusters
Vector rotor;
vector< double > mbefore, mafter; // fitting mismatches
vector< double > mstretch, mbend; //bond length change effect, residual (bend) effect
vector< Vector > mrot; // the rotors
vector< double > modrot; // size of vector
vector< double > indeg; // in degrees

bool readpol( string filename, vector< Cluster > &clus, vector< Cluster > &ghost ); // returns false if the read is bad.

int main(int argc, char *argv[]){
    
  string infilename;
  
  //cerr << "Args: " << argc << endl;
  if ( argc < 4 ) {
     stringstream whatn;
     whatn << argv[0];
     string myname = whatn.str();
     cerr << "Insufficient filenames given." << endl;
     cerr << "Please give two filenames with pol files for readin," << endl;
     cerr << " and a filename for the comparison output." << endl;
     cerr << " thus: " << myname << " in1 in2 out" << endl;
     cerr << " Will return rotations which fit the polyhedra of in2 (mobile) over those of in1 (target)," << endl;
     cerr << "  and the mismatches before and after fitting takes place." << endl;
     exit(0); 
  }
  
  stringstream what1, what2, what3;
  what1 << argv[1];
  what2 << argv[2];
  what3 << argv[3];
  
  string name1, name2, name3;
  name1 = what1.str();
  name2 = what2.str();
  name3 = what3.str();    

  //bool goodread = readpol ( infilename, clus1, ghost1 );
  bool good1 = readpol ( name1, clus1, ghost1 );
  bool good2 = readpol ( name2, clus2, ghost2 );
  
  if (!good1){
     cerr << "Problem with reading " << name1 << endl;
     exit(1);  
  }
  if (!good2){
     cerr << "Problem with reading " << name2 << endl;
     exit(1);  
  }        

  int n1 = clus1.size();
  int n2 = clus2.size();
  
  if ( n1 != n2 ){
       cerr << "Problem: disparate list sizes " << n1 << " , " << n2 << endl;
       exit(1);
  }

  //loop and solve cluster by cluster ...
  for ( int i = 0 ; i < n1 ; i++ ){
      int m1 = clus1.at(i).bond.size();
      int m2 = clus2.at(i).bond.size();
      
      if ( m1 != m2 ){
           cerr << "Problem: disparate cluster sizes at cluster " << i+1;
           cerr << " ; " << m1 << " , " << m2 << endl;
           exit(1);
      }
      
      //okay replace the stop with a warning flag and then check for reordering
      bool reorder = false;
      for (int j = 0 ; j < m1 ; j++ ){
          int a1 = clus1.at(i).members.at(j);
          int a2 = clus2.at(i).members.at(j);

          if ( a1 != a2 ) {
               cerr << "Problem: mismatched atom ids at cluster " << i+1;
               cerr << " , member " << j+1;
               cerr << " ; " << a1 << " , " << a2 << endl;
               reorder = true;
          }
          
      }
      
      //and now attempt reordering to match a1/a2 pairs
      //stop/fail only if the reordering is not possible
      //otherwise reorder the mobile (clus2) to match the target (clus1)
      if ( reorder ) {
         Cluster tempclus = clus2.at(i); // temporary copy to play with
         vector < bool > ticks ; //mark ids as found, final trap for any failures
         ticks.resize( m2 ); // one per done
         
         for ( int j = 0; j < m1 ; j++ ){
             ticks.at(j) = false; // set to true when good
         }
         for ( int j = 0; j < m1 ; j++ ){
             int a1 = clus1.at(i).members.at(j) ;
             for ( int k = 0 ; k < m2 ; k++ ){
                 int a2 = clus2.at(i).members.at(k);
                 if ( a1 == a2 ){
                    //match!
                    ticks.at(k) = true; // matched entry k2 to entry j1
                    tempclus.members.at(j) = a2; // id order now matches clus1.at(i)
                    tempclus.bond.at(j) = clus2.at(i).bond.at(k); // reorder bonds too
                    cerr << "Matched " << a1 << " to " << a2 << endl;
                 }   
             }
         }
         
         //did we succeed?
         for ( int j = 0; j < m1 ; j++ ){
             if ( ticks.at(j) ){
                  continue; // good, keep going
             }
             else{
                  cerr << "FATAL problem: no match for " << clus1.at(i).members.at(j) << endl;
                  exit(1);
             }
         }         
         clus2.at(i) = tempclus ; // overwrite with reordered values
      
      } // done reordering!
      
      //if we get here then the clusters should be a good match!
      // first do vector mismatch before rotation fitting;
      double mis = 0.0;
      Vector diff;
      
      for (int j = 0 ; j < m1 ; j++ ){
          diff = clus1.at(i).bond.at(j) - clus2.at(i).bond.at(j);
          mis += diff.sq();
          
      }      
      cerr << i+1 << " : Mismatch-square before fitting: " << mis << endl;
      mbefore.push_back( mis );

      rotor = nullvec;
      fitrotor( clus1.at(i), clus2.at(i), rotor );
      
      cerr << "Rotor " << rotor << endl;
      mrot.push_back( rotor );
      
      //some calculations:
      double rotsize = sqrt( rotor.sq() ); // magnitude
      modrot.push_back( rotsize );
      cerr << "Bivector size " << rotsize;
      
      // and convert to degrees
      // twice the sine of half the angle...
      double foo = rotsize * 0.5;
      double baz = asin( foo );
      foo = 2* baz * radtodeg;
      cerr << " -> " << foo << " degrees rotation angle." << endl;
      indeg.push_back( foo );
      
      //now mismatch after
      mis = 0.0;
      
      for (int j = 0 ; j < m1 ; j++ ){
          diff = clus1.at(i).bond.at(j) - clus2.at(i).bond.at(j);
          mis += diff.sq();
          
      }      
      cerr << i+1 << " : Mismatch-square after fitting: " << mis << endl; 
      mafter.push_back(mis);
      
      //and finally do the stretch versus bend breakdown.
      double stretches = 0.;
      for (int j = 0 ; j < m1 ; j++ ){
		  double l1 = sqrt( clus1.at(i).bond.at(j).sq() ); //bond length 1
		  double l2 = sqrt( clus2.at(i).bond.at(j).sq() ); //bond length 2
		  stretches += ( l1-l2)*(l1-l2); //square of change in bond length          
      }       
      
      cerr << i+1 << " : Stretching mismatch (bond-length-change component): " << stretches << endl;
      mstretch.push_back( stretches );
      double residual = mis - stretches;
      cerr << i+1 << " : Bending mismatch (non-bond-length-change component): " << residual << endl;      
      mbend.push_back( residual );
      
      
  }

  //write to output file
  ofstream outfile( name3.c_str() );
  outfile << "ID M2-before M2-after M2-bend M2-stretch Rx Ry Rz ModR InDegrees" << endl;
  for ( int i = 0 ; i < n1 ; i++ ){
      outfile << i+1 << " " << mbefore.at(i) << " ";
      outfile << mafter.at(i) << " ";
      outfile << mbend.at(i) << " ";
      outfile << mstretch.at(i) << " ";
      outfile <<mrot.at(i) << " ";
      outfile << modrot.at(i) << " ";
      outfile << indeg.at(i) << endl;
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
