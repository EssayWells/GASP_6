#include "basic_utils.h" // the essentials including random and twopi

#ifndef __VECTORS_H_INCLUDED__
#define __VECTORS_H_INCLUDED__


class Vector{
 public:
 double x,y,z;
 
 Vector(double x_in=0., double y_in=0., double z_in=0.){
  x = x_in; y = y_in; z = z_in;
 }
 
 Vector operator+ (const Vector& addme);

 Vector operator- (const Vector& addme);
 
  Vector operator* (const double& timesme);
 
   Vector operator/ (const double& timesme);
 
  friend ostream &operator<<(ostream &out, const Vector &c);
 
   double dot (const Vector &dotme);
 
   Vector cross (const Vector &crossme);
 
   Vector& operator+=(const Vector& plusme);
   Vector& operator-=(const Vector& lessme);
   Vector& operator*=(const double& timesme);
   Vector& operator/=(const double& divme);
    
   friend Vector operator*(double timesme, const Vector &avec);

   double sq(); //square of vector
   
   void norm();
};

Vector unit_vector( ); //random unit vector
Vector cubic_vector(); //random vector in a cube (-0.5 to +0.5 each component)

Vector wrapfracpos( Vector fracin );
Vector wrapfracdel( Vector fracin );

const Vector nullvec=Vector(0.,0.,0.);

Vector first3( vector< double > &params );
Vector second3( vector< double > &params );
vector< double > sixfromtwovecs( Vector &sidein, Vector &anglein );
#endif
