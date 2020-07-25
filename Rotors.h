#include "basic_utils.h"
#include "Vectors.h"

#ifndef __ROTORS_H_INCLUDED_
#define __ROTORS_H_INCLUDED_

double Xval(Vector &b);
double epsx(Vector &pq, double pqdash_x, Vector &b, double myX);
double epsy(Vector &pq, double pqdash_y, Vector &b, double myX);
double epsz(Vector &pq, double pqdash_z, Vector &b, double myX);
double depsxdbx(Vector &pq, Vector &b, double myX);
double depsxdby(Vector &pq, Vector &b, double myX);
double depsxdbz(Vector &pq, Vector &b, double myX);
double depsydbx(Vector &pq, Vector &b, double myX);
double depsydby(Vector &pq, Vector &b, double myX);
double depsydbz(Vector &pq, Vector &b, double myX);
double depszdbx(Vector &pq, Vector &b, double myX);
double depszdby(Vector &pq, Vector &b, double myX);
double depszdbz(Vector &pq, Vector &b, double myX);
double deps2dbx(Vector &pq, Vector &pqdash, Vector &b);
double deps2dby(Vector &pq, Vector &pqdash, Vector &b);
double deps2dbz(Vector &pq, Vector &pqdash, Vector &b);

struct Cluster{
  public:
         
         vector< int > members;
         Vector cpos;
         vector< Vector > bond; // all bond vectors here!
	 Cluster(){}
         
         void randomorient(); // randomly reorient to a new basis set
};

bool isincluster( int id , Cluster &clust );
void fitrotor( Cluster &target, Cluster &mobile, Vector &rotor);
Vector rotated( Vector rotateme, Vector rotor, double myX );
void rotate( Cluster &target, Vector rotor);
Cluster uniclust( Cluster &a, Cluster &b );
const Cluster nullcluster; // an empty cluster object for placeholding


struct Poly{
  public:

         string c_species;
         string v_species;
         int shape;
         double bondlength;
         vector< int > members;
         Vector cpos;
         vector< Vector > bond; // all bond vectors here! 
         int type;
         Poly(){
			 c_species=""; v_species=""; shape=0; bondlength=0.0; type=0;
		 }
         
         
};

const Poly nullpoly; // empty poly for placeholding
Cluster clusterfrompoly( Poly &poly );
void putclusterinpoly( Cluster &clust, Poly &poly);
void perfectghost( Poly &aghost ); // after poly setup, do the geometry right
void buildghost( Poly &aghost ); // create geometry from scratch

#endif

