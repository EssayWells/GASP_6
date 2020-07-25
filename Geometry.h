#include "basic_utils.h"
#include "Vectors.h"
#include "Rotors.h"


#ifndef __GEOMETRY_H_INCLUDED_
#define __GEOMETRY_H_INCLUDED_

struct Bondspec{
  public:
         
         string speciesA;
         vector< string > speciesB;
         double within;      
         Bondspec(){
			 speciesA=""; within=0.;
		 }   
};

const Bondspec nullbondspec;
//vector< Bondspec > bondspec; // when we get them


struct Polyspec{
  public:

         string c_species;
         string v_species;
         int shape;
         double bondlength;
         int nmembers; // number of this type of poly; needed?
         Polyspec(){
			 c_species=""; v_species=""; shape=0; bondlength=0.1; nmembers=0;
		 }
};

const Polyspec nullpolyspec;
//vector< Polyspec > polyspec;


class Bondline{
  public:
         
         int id1;
         int id2;
         int bars;
         string ele1;
         string ele2;
         Bondline(int id1in=0, int id2in=0, int barsin=5, string ele1in="", string ele2in=""){
			 id1=id1in; id2=id2in; bars=barsin; ele1=ele1in; ele2=ele2in;
		 }
		 string printbond();
         
};


struct Contact{
       public:
              Vector dr; // relative position of contacting atom
              double l; // contact distance
              int id1; // central atom
              int id2; // other atom
              Contact(){
				  l=0.1; id1=0; id2=0;
			  }
};
  
void writebonding( string filename, vector< Bondline > &bl );
bool readbonds( string filename, vector< Bondline > &bl );



#endif
