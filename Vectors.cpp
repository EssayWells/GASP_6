#include "Vectors.h"


 Vector Vector::operator+ (const Vector& addme){
        Vector temp;
        temp.x=x+addme.x;
        temp.y=y+addme.y; 
        temp.z=z+addme.z;
        return(temp);
 }

 Vector Vector::operator- (const Vector& addme){
        Vector temp;
        temp.x=x-addme.x;
        temp.y=y-addme.y; 
        temp.z=z-addme.z;
        return(temp);
 }
 
  Vector Vector::operator* (const double& timesme){
        Vector temp;
        temp.x=x*timesme;
        temp.y=y*timesme; 
        temp.z=z*timesme;
        return(temp);
 }
 
   Vector Vector::operator/ (const double& timesme){
        Vector temp;
	double nome = 1.0/timesme;
        temp.x=x*nome;
        temp.y=y*nome; 
        temp.z=z*nome;
        return(temp);
 }
 
ostream &operator<<(ostream &out, const Vector &c){
        out << fixed << showpoint << setw(12) << setprecision(8) << c.x;
        out << " ";
        out << fixed << showpoint << setw(12) << setprecision(8) << c.y;
        out << " ";
        out << fixed << showpoint << setw(12) << setprecision(8) << c.z;
        return out;
 }
 
   double Vector::dot (const Vector &dotme){
          return x*dotme.x  + y*dotme.y + z*dotme.z;
   }
 
   Vector Vector::cross (const Vector &crossme){
          Vector temp;
          temp.x = y*crossme.z  - z*crossme.y;
          temp.y = z*crossme.x - x*crossme.z;
          temp.z = x*crossme.y - y*crossme.x;
          return temp;
   }
    

   double Vector::sq(){
          return x*x + y*y + z*z; // vector-squared shortcut!
   }
   
   void Vector::norm(){
        double d2 = x*x + y*y + z*z;
        if ( d2 < 1.0e-100 ) return; // do not norm a null vector!
        double normer = 1.0/sqrt( d2 ) ;
        x *= normer;
        y *= normer;
        z *= normer;
        return; // the vector has now been normalised
   }

 //not sure about Vector:: for this one, check.
Vector operator*(double timesme, const Vector &avec){
 Vector result = avec;
 result *= timesme;
 return result;            
} // this overloads double * Vector
//we used a member function to overload Vector * double

Vector& Vector::operator+=(const Vector& plusme)
{
        x += plusme.x;
        y += plusme.y;
        z += plusme.z;
    return *this;
};

 
Vector& Vector::operator-=(const Vector& lessme)
{
        x -= lessme.x;
        y -= lessme.y;
        z -= lessme.z;
    return *this;
};

Vector& Vector::operator*=(const double& timesme){
        x *= timesme;
        y *= timesme;
        z *= timesme;        
        return *this;
};

Vector& Vector::operator/=(const double& divme){
        x /= divme;
        y /= divme;
        z /= divme;
        return *this;
};


Vector unit_vector( ){
  double myphi, mycostheta, mysintheta, range;
  Vector res(0,0,0);
      range = mytwopi;
      myphi = range*getRandom(); // 0 to 2pi
      range = 2;
      mycostheta = ( range*getRandom() ) - 1.0; //-1 to 1
      mysintheta = sqrt(1 - pow(mycostheta,2));

      res.x = mysintheta * cos(myphi);
      res.y = mysintheta * sin(myphi); 
      res.z = mycostheta;
  
  return res;
} 

Vector cubic_vector(){
	double rx, ry, rz;
	rx= getRandom() - 0.5;
	ry= getRandom() - 0.5;	
	rz= getRandom() - 0.5;
	Vector result = Vector( rx, ry, rz );
	return result;
}

//
// wrapfracpos ensures that a fractional coordinate is not less than zero, nor equal to or greater than 1.
//
Vector wrapfracpos( Vector fracin ){
       Vector result;
       result = fracin;
       while ( result.x >= 1.0 ) result.x -= 1.0;
       while ( result.x <  0.0 ) result.x += 1.0;
       while ( result.y >= 1.0 ) result.y -= 1.0;
       while ( result.y <  0.0 ) result.y += 1.0;
       while ( result.z >= 1.0 ) result.z -= 1.0;
       while ( result.z <  0.0 ) result.z += 1.0;
       return result;       
}

//
// wrapfracdel ensures that a fractional difference does not cross more than half a cell.
//
Vector wrapfracdel( Vector fracin ){
       Vector result;
       result = fracin;
       while ( result.x >= 0.5 ) result.x -= 1.0;
       while ( result.x <  -0.5 ) result.x += 1.0;
       while ( result.y >= 0.5 ) result.y -= 1.0;
       while ( result.y <  -0.5 ) result.y += 1.0;
       while ( result.z >= 0.5 ) result.z -= 1.0;
       while ( result.z <  -0.5 ) result.z += 1.0;
       return result;       
}


Vector first3( vector< double > &params ){
	if ( params.size() < 3 ){
		return Vector( 0.,0.,0.);
	}
	Vector result = Vector( params.at(0), params.at(1), params.at(2) );
	return result;
}
Vector second3( vector< double > &params ){
		if ( params.size() < 6 ){
		return Vector( 0.,0.,0.);
	}
	Vector result = Vector( params.at(3), params.at(4), params.at(5) );
	return result;
}
vector< double > sixfromtwovecs( Vector &sidein, Vector &anglein ){
	vector< double > result;
	result.resize(6);
	result.at(0) = sidein.x;
	result.at(1) = sidein.y;
	result.at(2) = sidein.z;
	result.at(3) = anglein.x;
	result.at(4) = anglein.y;
	result.at(5) = anglein.z;
	return result;
}


