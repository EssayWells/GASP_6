#include "Rotors.h" //which brings basic utils and Vectors


//////////////////////////////////////////////////////////////////////////////////
// Decription: The geometric algebra equations for rotating bonds.
////////////////////////////////////////////////////////////////////////////////
double Xval(Vector &b){
  return sqrt( 1.0 - 0.25 * (b.x*b.x + b.z*b.z + b.y*b.y)); 
}

////////////////////////////////////////////////////////////////////////////////
double epsx(Vector &pq, double pqdash_x, Vector &b, double myX){
  double eps_x;
  
  eps_x = pq.x*(1.0 - 0.5*(b.y*b.y+ b.z*b.z))
    + pq.y*(-myX*b.z + 0.5*b.x*b.y)
    + pq.z*(myX*b.y + 0.5*b.x*b.z)
    - pqdash_x;
  
  return eps_x;
}

////////////////////////////////////////////////////////////////////////////////
double epsy(Vector &pq, double pqdash_y, Vector &b, double myX){
  double eps_y;
  
  eps_y = pq.y*(1.0 - 0.5*(b.x*b.x + b.z*b.z))
    + pq.z*(-myX*b.x + 0.5*b.y*b.z)
    + pq.x*(myX*b.z + 0.5*b.x*b.y)
    - pqdash_y;
  
  return eps_y;
}

////////////////////////////////////////////////////////////////////////////////
double epsz(Vector &pq, double pqdash_z, Vector &b, double myX){
  double eps_z;
  eps_z = pq.z*(1.0 - 0.5*(b.x*b.x + b.y*b.y))
    + pq.x*(-myX*b.y + 0.5*b.x*b.z)
    + pq.y*(myX*b.x + 0.5*b.y*b.z)
    - pqdash_z;
  
  return eps_z;
}

/* deps{m}db{n}(pq,b) is the differential of eps{m} by b{n} */

////////////////////////////////////////////////////////////////////////////////
double depsxdbx(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.y*( (0.25*b.x*b.z/myX) + 0.5*b.y)
  //+ pq.z*( (-0.25*b.x*b.y/myX) + 0.5*b.z);

  double res = pq.y*( (0.25*b.x*b.z/myX) + 0.5*b.y) 
    + pq.z*( (-0.25*b.x*b.y/myX) + 0.5*b.z);
  
  return res;
}

////////////////////////////////////////////////////////////////////////////////
double depsxdby(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*(-b.y)
  //+ pq.y*( (0.25*b.y*b.z/myX) + 0.5*b.x)
  //+ pq.z*(myX - (0.25*b.y*b.y/myX));
  
  double res = pq.x*(-b.y)
    + pq.y*( (0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.z*(myX - (0.25*b.y*b.y/myX));
  
  return res;
}

////////////////////////////////////////////////////////////////////////////////
double depsxdbz(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*(-b.z)
  //+ pq.z*( (-0.25*b.y*b.z/myX) + 0.5*b.x)
  //+ pq.y*(-myX + (0.25*b.z*b.z/myX));
  
  double res = pq.x*(-b.z)
    + pq.z*( (-0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.y*(-myX + (0.25*b.z*b.z/myX));

  return res;
}

////////////////////////////////////////////////////////////////////////////////
double depsydbx(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*((-0.25*b.x*b.z/myX) + 0.5*b.y)
  //+ pq.y *(-b.x)
  //+ pq.z*(-myX + (0.25*b.x*b.x/myX));
  
  double res = pq.x*((-0.25*b.x*b.z/myX) + 0.5*b.y)
    + pq.y *(-b.x)
    + pq.z*(-myX + (0.25*b.x*b.x/myX)); 
  return res;
}

////////////////////////////////////////////////////////////////////////////////
double depsydby(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*((-0.25*b.y*b.z/myX) + 0.5*b.x)
  //+ pq.z*((0.25*b.x*b.y/myX) + 0.5*b.z);

  double res = pq.x*((-0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.z*((0.25*b.x*b.y/myX) + 0.5*b.z);
  return res;
}

////////////////////////////////////////////////////////////////////////////////
double depsydbz(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*(myX - (0.25*b.z*b.z/myX))
  //+ pq.y*(-b.z)
  //+ pq.z*((0.25*b.x*b.z/myX) + 0.5*b.y);

  double res = pq.x*(myX - (0.25*b.z*b.z/myX))
    + pq.y*(-b.z)
    + pq.z*((0.25*b.x*b.z/myX) + 0.5*b.y);
  return res;
}

////////////////////////////////////////////////////////////////////////////////
double depszdbx(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*((0.25*b.x*b.y/myX) + 0.5*b.z)
  //+ pq.y*(myX - (0.25*b.x*b.x/myX))
  //+ pq.z*(-b.x);
  
  double res = pq.x*((0.25*b.x*b.y/myX) + 0.5*b.z)
    + pq.y*(myX - (0.25*b.x*b.x/myX))
    + pq.z*(-b.x);
  return res;
} 

////////////////////////////////////////////////////////////////////////////////
double depszdby(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*(-myX + (0.25*b.y*b.y/myX))
  //+ pq.y*((-0.25*b.x*b.y/myX) + 0.5*b.z)
  //+ pq.z*(-b.y);
  
  double res = pq.x*(-myX + (0.25*b.y*b.y/myX))
    + pq.y*((-0.25*b.x*b.y/myX) + 0.5*b.z)
    + pq.z*(-b.y);
  return res;
} 

////////////////////////////////////////////////////////////////////////////////
double depszdbz(Vector &pq, Vector &b, double myX){
  //double res;
  //res = pq.x*((0.25*b.y*b.z/myX) + 0.5*b.x)
  //+ pq.y*((-0.25*b.x*b.z/myX) + 0.5*b.y);
  
  double res = pq.x*((0.25*b.y*b.z/myX) + 0.5*b.x)
    + pq.y*((-0.25*b.x*b.z/myX) + 0.5*b.y); 
  return res;
} 

/*The functions deps2dbn are the mismatch */
/* eps squared differentiated by rotor component bn */
/* These are used to minimise eps squared wrt the rotor b */

////////////////////////////////////////////////////////////////////////////////
double deps2dbx(Vector &pq, Vector &pqdash, Vector &b){
  double myX, res;
  myX = Xval(b);
  
  res = 2*epsx(pq, pqdash.x, b, myX)*depsxdbx(pq, b, myX)
    + 2*epsy(pq, pqdash.y, b, myX)*depsydbx(pq, b, myX)
    + 2*epsz(pq, pqdash.z, b, myX)*depszdbx(pq, b, myX);
  
  return res;
}

////////////////////////////////////////////////////////////////////////////////
double deps2dby(Vector &pq, Vector &pqdash, Vector &b){
  double myX, res;
  myX = Xval(b);
  
  res = 2*epsx(pq, pqdash.x, b, myX)*depsxdby(pq, b, myX)
    + 2*epsy(pq, pqdash.y, b, myX)*depsydby(pq, b, myX)
    + 2*epsz(pq, pqdash.z, b, myX)*depszdby(pq, b, myX);
  
  return res;
}

////////////////////////////////////////////////////////////////////////////////
double deps2dbz(Vector &pq, Vector &pqdash, Vector &b){
  double myX, res;
  myX = Xval(b);
  
  res = 2*epsx(pq, pqdash.x, b, myX)*depsxdbz(pq, b, myX)
      + 2*epsy(pq, pqdash.y, b, myX)*depsydbz(pq, b, myX)
      + 2*epsz(pq, pqdash.z, b, myX)*depszdbz(pq, b, myX);
  
  return res;
}

//header includes Cluster data type
void Cluster::randomorient(){
	//first create a new basis set of unit vectors
	Vector basis1 = unit_vector(); 
	Vector basis3 = unit_vector(); //i, k' random
	Vector basis2 = basis3.cross( basis1 ); // create RH j
	basis2.norm(); // make it unit
	basis3 = basis1.cross( basis2 ); // recreate RH k as unit vector by construction
	//then multiply bonds by the new basis
	for ( int i = 0 ; i < bond.size(); i++ ){
		Vector newb = nullvec; // zero at start
		newb = basis1 * bond.at(i).x  + basis2 * bond.at(i).y + basis3 * bond.at(i).z;
		bond.at(i) = newb;
	}
	//and we are done
}


//
//isincluster returns true if id is a member of cluster
//
bool isincluster( int id , Cluster &clust ){
     for ( int i = 0 ; i < clust.members.size() ; i++ ){
         if ( id == clust.members.at(i) ) return true;
     }
     return false;
}

void fitrotor( Cluster &target, Cluster &mobile, Vector &rotor){
     double step;
     double smallgsq = 1e-20; // small gradient squared = fitted. Tune?
     double largersq = 0.01; // big rotor = time to update basis during fit
     double dot;
     Vector rot ; // use this internally
     //relevant vectors in in target and mobile
     
     Vector g0, g1; // gradients of mismatches
     double g02, g12; // grad-squared
     double alpha; // used in secant
     Vector trial; // trial vector = needed?
     Vector tempvec; // working space
     
     bool update_base = false;
     int maxc=50; // shold have fitted well within 50 cycles!
     int msize = mobile.bond.size();
     
     Vector diffvec;
     double differ = 0.0; // just checking - mismatch before/after fitting
     
     rot = nullvec; //should I allow rot=input rotor
     step = 1.0/( 10.0 * msize * msize); // a rough stepsize in rotor space
     
     for ( int k = 0 ; k < maxc; k++ ){
         //cerr << "Fitting cycle " << k << endl;
         bool leave = false; // use to flag fitting end
         g0 = nullvec;
         g1 = nullvec;
         for ( int p = 0 ; p < msize ; p++){
             g0.x +=deps2dbx( mobile.bond.at(p), target.bond.at(p), rot);
             g0.y +=deps2dby( mobile.bond.at(p), target.bond.at(p), rot);
             g0.z +=deps2dbz( mobile.bond.at(p), target.bond.at(p), rot);
         }
         g02 = g0.sq();
         //cerr << "Rotor = " << rot << "; ";
         //cerr << "Gradsquared = " << g02 << endl;
         if ( g02 < smallgsq ) break; // well fitted!
         
         tempvec = g0 * step; //step along gradient
         dot = tempvec.sq();
         if ( dot > largersq ){
              //trap large rotor for stability
              tempvec *= ( largersq/dot ); // scale down
         }
         trial = rot - tempvec; // trial rotor value
         for ( int p = 0 ; p < msize ; p++){
             g1.x +=deps2dbx( mobile.bond.at(p), target.bond.at(p), trial);
             g1.y +=deps2dby( mobile.bond.at(p), target.bond.at(p), trial);
             g1.z +=deps2dbz( mobile.bond.at(p), target.bond.at(p), trial);
         }
         g12 = g1.sq();
         // g1 is gradient when rotor = trial
         if ( g12 < smallgsq){
              leave = true;
              rot = trial;
              //cerr << "Rotor = " << rot << "; ";
              //cerr << "Gradsquared = " << g12 << endl;
         }
         if ( leave ) break; // exit loop, with rotor at rot 
         
         //use secant if appropriate
         if ( k > 1 && g12 < g02 ){
              alpha = 1.0 / ( 1.0 - ( g0.dot(g1)/g02 ) );
              tempvec *= alpha; // better step estimation
         }
         
         //update rot value:
         dot = tempvec.sq();
         if ( dot > largersq){
              tempvec *= largersq/dot; // scale down big steps
         }
         rot -= tempvec; // updated rot
         dot = rot.sq();
         if ( dot > largersq )update_base = true; //start rotating ghost to avoid instability
         
         //experiment:
                      //update_base = true;
         if ( update_base){
            rotate( mobile, rot);
            rotor += rot; //running update on "rotor" argument!
            rot = nullvec;
            update_base = false;
         }
         
     }
     
     rotate( mobile, rot); // rotate by final rotor
     rotor += rot; // update returned rotor
          
     return;
}

//
//rotated takes a vector and a rotor (plus Xval) and rotates the vector
//
Vector rotated( Vector rotateme, Vector rotor, double myX ){
       Vector result;
          result.x = epsx(rotateme, 0.0, rotor, myX); 
          result.y = epsy(rotateme, 0.0, rotor, myX);
          result.z = epsz(rotateme, 0.0, rotor, myX);          
       return result;
} // END OF ROTATED

//
//rotate applies a rotor to a cluster, using "rotated" on each bond
//
void rotate( Cluster &target, Vector rotor){
     Vector tempvec;
     double anX = Xval( rotor );
     for ( int i = 0 ; i < target.bond.size(); i++ ){
         tempvec = rotated( target.bond.at(i), rotor, anX ) ;
         target.bond.at(i) = tempvec;    
     } 
     return;
} // END OF ROTATE

//
//union unites the membership of a and b
//
Cluster uniclust( Cluster &a, Cluster &b ){
        Cluster result;
        result = a ;
        int bsize = b.members.size();
        for ( int i = 0; i < bsize; i++){
            int j = b.members.at(i);
            if ( ! isin( j, result.members ) ){
                 result.members.push_back( j );
                 result.bond.push_back( b.bond.at(i) ); //just a placeholder
            }
        }
        return result;
}



//
//puts poly geometry into a Cluster object ready for fitting
//
Cluster clusterfrompoly( Poly &poly ){
        Cluster result;
        result.members = poly.members;
        result.cpos = poly.cpos;
        result.bond = poly.bond;
        return result;
}

//
//puts the cluster geometry into the given Poly
//
void putclusterinpoly( Cluster &clust, Poly &poly){
	 poly.members = clust.members;
     poly.cpos = clust.cpos;
     poly.bond = clust.bond;
     return;
}

//
//perfectghost reshapes a ghost poly into the ideal geometry
//
void perfectghost( Poly &aghost ){
     //aghost.shape
     //aghost.bondlength
     //keep cpos
     //adjust bond array
     Vector basis1, basis2, basis3, cross, swapvec;
     Vector tempvec, workvec; //working space if needed
     double test, most, least, size, dot;
     bool rh; // right handed basis set?
     int up, down, left, right, towards, away; // for arranging octahedral bonds
     
     double one3 = 1.0/3.0;
     double root2 = sqrt(2.0);
     double root3 = sqrt(3.0);
     double oneroot3 = 1.0/root3;
     
     rh=false;
     test=most=least=size=dot=0;
     up=down=left=right=towards=away=0;
     
     if ( aghost.shape == 1 ){
          //tetrahedron case
          //centre and zero:
          workvec = aghost.bond.at(0);
          for ( int j = 0 ; j < 5 ; j++ ){
              aghost.bond.at(j) -= workvec;
          }
          basis1 = aghost.bond.at(1);
          basis2 = aghost.bond.at(2);
          basis3 = aghost.bond.at(3);
          
          cross = basis1.cross( basis2 );
          test = cross.dot( basis3 );
          if ( test > 0.0 ) rh = true; // handedness
          basis1.norm();
          basis3 = cross;
          basis3.norm();
          basis2 = basis3.cross( basis1 ); // orthonormal
          
          aghost.bond.at(1) = basis1;
          workvec = basis1*( - one3 ) + basis2 * ( 2.0 * root2 * one3 ) ;
          aghost.bond.at(2) = workvec;
          //build rh bond first
          workvec = basis1* ( - one3 ) - basis2 *root2 * one3 + basis3 * root2 * oneroot3;
          aghost.bond.at(3) = workvec;
          workvec = basis1* ( - one3 ) - basis2 *root2 * one3 - basis3 * root2 * oneroot3;
          aghost.bond.at(4) = workvec;
          //swap if not rh
          if ( !rh ){
               swapvec = aghost.bond.at(3);
               aghost.bond.at(3) = aghost.bond.at(4);
               aghost.bond.at(4) = swapvec;
          }
          //expand to right size
          for ( int j = 1 ; j < 5 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }          
          //done tetrahedron geometry!
     }
     if ( aghost.shape == 2 ){
          //octahedron case
          //centre and zero:
          workvec = aghost.bond.at(0);
          for ( int j = 0 ; j < 7 ; j++ ){
              aghost.bond.at(j) -= workvec;
          }          
          
          up = 1; // first bond defines up direction
          workvec = aghost.bond.at(1);
          down = 0; // for now!
          least = test = 0.0;
          for ( int j = 2 ; j < 7 ; j++ ){
              test = aghost.bond.at(j).dot( workvec );
              if ( test < least ){
                   down = j; least = test;
              }
          }           
          //now up = 1, down = j
          
          if ( down == 2 ){
               left = 3;
          }
          else{
               left = 2;
          }
          //defined left direction
          least = test = 0.0;
          workvec = aghost.bond.at(left);
          for ( int j = 3 ; j < 7 ; j++ ){
              if ( j == down ) continue;
              if ( j == left ) continue;
              test = aghost.bond.at(j).dot( workvec );
              if ( test < least ){
                   right = j; least = test;
              }
          }           
          // now got up down left and right!
          
          cross = aghost.bond.at(right).cross( aghost.bond.at( up ) );
          test = 0.0;
          workvec = cross;
          for ( int j = 3 ; j < 7 ; j++ ){
              if ( j == down ) continue;
              if ( j == left ) continue;
              if ( j == right ) continue;
              test = aghost.bond.at(j).dot( workvec );
              if ( test > 0.0 ){
                   towards = j;
              }
              else {
                   away = j;
              }
          }           
          //fully sorted octahedron!
          
          basis1 = aghost.bond.at(up);
          basis1.norm();
          aghost.bond.at(up) = basis1;
          aghost.bond.at(down) = basis1 * -1.0;
          
          basis2 = aghost.bond.at( right );
          dot = basis1.dot( basis2 );
          basis2 -= basis1*dot;
          basis2.norm();
          aghost.bond.at( right ) = basis2;
          aghost.bond.at( left ) = basis2 * -1.0;
          
          cross = basis1.cross( basis2 );
          aghost.bond.at( away ) = cross;
          aghost.bond.at( towards ) = cross *-1.0;
          
          //expand to right size
          for ( int j = 1 ; j < 7 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }           
     }     
     if ( aghost.shape == 3 ){
          //triangle case
          //centre and zero:
          workvec = aghost.bond.at(0);
          for ( int j = 0 ; j < 4 ; j++ ){
              aghost.bond.at(j) -= workvec;
          }
          basis1 = aghost.bond.at(1);
          basis2 = aghost.bond.at(2);
          cross = basis2.cross( basis1 );
          basis1.norm();
          basis3 = cross;
          basis3.norm();
          basis2 = basis1.cross( basis3 ); // orthonormal
          
          aghost.bond.at(1) = basis1;
          workvec = basis1*( - 0.5 ) + basis2 * ( root3 * 0.5 ) ;
          aghost.bond.at(2) = workvec;
          //build rh bond first
          workvec = basis1*( - 0.5 ) - basis2 * ( root3 * 0.5 ) ;
          aghost.bond.at(3) = workvec;

          //expand to right size
          for ( int j = 1 ; j < 4 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }          
          //done triangle geometry!
     }     
     if ( aghost.shape == 4 ){
          //square case
          //centre and zero:
          workvec = aghost.bond.at(0);
          for ( int j = 0 ; j < 5 ; j++ ){
              aghost.bond.at(j) -= workvec;
          }          
          
          up = 1; // first bond defines up direction
          workvec = aghost.bond.at(1);
          down = 0; // for now!
          least = test = 0.0;
          for ( int j = 2 ; j < 5 ; j++ ){
              test = aghost.bond.at(j).dot( workvec );
              if ( test < least ){
                   down = j; least = test;
              }
          }           
          //now up = 1, down = j
          
          if ( down == 2 ){
               left = 3;
               right = 4;
          }
          else if ( down == 3 ){
               left = 2;
               right = 4;
          }
          else {
               //down == 4
               left = 2;
               right = 3;               
          }
          //defined left direction
          // now got up down left and right!
          //fully sorted square!
          
          basis1 = aghost.bond.at(up);
          basis1.norm();
          aghost.bond.at(up) = basis1;
          aghost.bond.at(down) = basis1 * -1.0;
          
          basis2 = aghost.bond.at( right );
          dot = basis1.dot( basis2 );
          basis2 -= basis1*dot; // made perpendicular
          basis2.norm();
          aghost.bond.at( right ) = basis2;
          aghost.bond.at( left ) = basis2 * -1.0;
          //expand to right size
          for ( int j = 1 ; j < 5 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }           
     }     

     if ( aghost.shape == 5 ){
          //bar case
          //easy case: bonds are equal and opposite by construction!
          for ( int j = 0 ; j < 2 ; j++ ){
              aghost.bond.at(j).norm();
              aghost.bond.at(j) *= 0.5*aghost.bondlength;
          }
     }
     
     return;
     //test report here:
     for ( int j = 0 ; j < aghost.bond.size() ; j++ ) {
         size = sqrt ( aghost.bond.at(j).sq() );
         cerr << "Ghost vector: " << aghost.bond.at(j) << " , (" << size << ")" << endl;
     }
            
            
     return;
}


//
//buildghost creates the ideal geometry
//
void buildghost( Poly &aghost ){
     //aghost.shape
     //aghost.bondlength
     Vector basis1, basis2, basis3;
     Vector tempvec, workvec; //working space if needed
     
     double one3 = 1.0/3.0;
     double root2 = sqrt(2.0);
     double root3 = sqrt(3.0);
     double oneroot3 = 1.0/root3;
     
     //create a local basis set
     basis1 = unit_vector();
     basis3 = unit_vector(); // random start
     basis2 = basis3.cross( basis1 );
     basis2.norm();
     basis3 = basis1.cross( basis2 ); 
     //built local rh basis set with random orientation
     
     if ( aghost.shape == 1 ){
          //tetrahedron case
          aghost.members.resize(5);
          aghost.bond.resize(5);
          aghost.bond.at(0) = nullvec; // just setting centre
          for ( int j = 0 ; j < 5 ; j++ ){
              aghost.members.at(j) = j; // dummy 0-4
          }
          aghost.bond.at(1) = basis1;
          workvec = basis1*( - one3 ) + basis2 * ( 2.0 * root2 * one3 ) ;
          aghost.bond.at(2) = workvec;
          //build rh bond first
          workvec = basis1* ( - one3 ) - basis2 *root2 * one3 + basis3 * root2 * oneroot3;
          aghost.bond.at(3) = workvec;
          workvec = basis1* ( - one3 ) - basis2 *root2 * one3 - basis3 * root2 * oneroot3;
          aghost.bond.at(4) = workvec;
          //expand to right size
          for ( int j = 1 ; j < 5 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }          
          //done tetrahedron geometry!
     }
     if ( aghost.shape == 2 ){
          //octahedron case
          aghost.members.resize(7);
          aghost.bond.resize(7);
          aghost.bond.at(0) = nullvec; // just setting centre
          for ( int j = 0 ; j < 7 ; j++ ){
              aghost.members.at(j) = j; // dummy 0-6
          }         
          
		  //set up one bond along each basis vector
		  aghost.bond.at(1) = basis1;
		  aghost.bond.at(2) = basis2;
		  aghost.bond.at(3) = basis3;
		  aghost.bond.at(4) = -1 * basis1;
		  aghost.bond.at(5) = -1 * basis2;
		  aghost.bond.at(6) = -1 * basis3;
          //expand to right size
          for ( int j = 1 ; j < 7 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }           
     }     
     if ( aghost.shape == 3 ){
          //triangle case
          aghost.members.resize(4);
          aghost.bond.resize(4);
          aghost.bond.at(0) = nullvec; // just setting centre
          for ( int j = 0 ; j < 4 ; j++ ){
              aghost.members.at(j) = j; // dummy 0-3
          } 
          
          aghost.bond.at(1) = basis1;
          workvec = basis1*( - 0.5 ) + basis2 * ( root3 * 0.5 ) ;
          aghost.bond.at(2) = workvec;
          workvec = basis1*( - 0.5 ) - basis2 * ( root3 * 0.5 ) ;
          aghost.bond.at(3) = workvec;

          //expand to right size
          for ( int j = 1 ; j < 4 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }          
          //done triangle geometry!
     }     
     if ( aghost.shape == 4 ){
          //square case
          aghost.members.resize(5);
          aghost.bond.resize(5);
          aghost.bond.at(0) = nullvec; // just setting centre
          for ( int j = 0 ; j < 5 ; j++ ){
              aghost.members.at(j) = j; // dummy 0-4
          }          
          
          //a square in the xy plane
		  aghost.bond.at(1) = basis1;
		  aghost.bond.at(2) = basis2;
		  aghost.bond.at(3) = -1 * basis1;
		  aghost.bond.at(4) = -1 * basis2;          

          //expand to right size
          for ( int j = 1 ; j < 5 ; j++ ){
              aghost.bond.at(j) *= aghost.bondlength;
          }           
     }     

     if ( aghost.shape == 5 ){
          //bar case
          aghost.members.resize(2);
          aghost.bond.resize(2);
          aghost.bond.at(0) = nullvec; // just setting centre
          for ( int j = 0 ; j < 2 ; j++ ){
              aghost.members.at(j) = j; // dummy 0-1
          }          
          
          //easy case: bonds are equal and opposite by construction         
		  aghost.bond.at(0) = basis1;
		  aghost.bond.at(1) = -1 * basis1;

          for ( int j = 0 ; j < 2 ; j++ ){
              aghost.bond.at(j) *= 0.5*aghost.bondlength;
          }
     }
     
     return;
}



