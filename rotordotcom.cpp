#include "Vectors.h"

vector< Vector > rot1, rot2 ; // the rotors in each .com file
double size1, size2, scalarprod ; // some magnitudes
bool readcom( string filename, vector< Vector > &rot); // returns false if the read is bad.
double dotme ( vector< Vector > &r1 , vector< Vector > &r2 , double &mod1, double &mod2 ); // scalar product
double dotdot (vector< Vector > &r1 ); // returns sum of components squared, is all.

int main(int argc, char *argv[]){
    
  string infilename;
  
  //cerr << "Args: " << argc << endl;
  if ( argc < 3 ) {
     stringstream whatn;
     whatn << argv[0];
     string myname = whatn.str();
     cerr << "Insufficient filenames given." << endl;
     cerr << "Please give two filenames with com files for readin" << endl;
     cerr << " thus: " << myname << " in1 in2" << endl;
     cerr << "The code will find the scalar product of the rotors in the two files." << endl;
     cerr << "This defines the similarity of two folding mechanisms" << endl;
     cerr << "Results written to standard output " << endl;
     cerr << "Use redirect > to put it in a file." << endl;
     exit(0); 
  }
    
  stringstream what1, what2;
  what1 << argv[1];
  what2 << argv[2];
  
  string name1, name2;
  name1 = what1.str();
  name2 = what2.str();  

  bool good1 = readcom ( name1, rot1 );
  bool good2 = readcom ( name2, rot2 );
  
  if (!good1){
     cerr << "Problem with reading " << name1 << endl;
     exit(1);  
  }
  if (!good2){
     cerr << "Problem with reading " << name2 << endl;
     exit(1);  
  }        

  int n1 = rot1.size();
  int n2 = rot2.size();
  
  if ( n1 != n2 ){
       cerr << "Problem: disparate list sizes " << n1 << " , " << n2 << endl;
       exit(1);
  }

  scalarprod = dotme ( rot1 , rot2 , size1 , size2 );
  
  cout << "Rotors 1: " << n1 << " rotors, magnitude " << size1 << endl;
  cout << "Rotors 2: " << n2 << " rotors, magnitude " << size2 << endl; 
  cout << "Scalar product: " << scalarprod << endl;
  
  double tiny = 1e-20; // piffling
  if ( size1 < tiny || size2 < tiny ) {
       cout << "Pseudo cos(theta): 0.0" << endl; // no angle to define for zero input vector
  }
  else{
       double ct = scalarprod / ( size1 * size2 ) ;
       cout << "Pseudo cos(theta): " << ct << endl;     
  }  
        
  exit(0); //final exit
} // END OF MAIN



double dotme ( vector< Vector > &r1 , vector< Vector > &r2 , double &mod1, double &mod2 ){
       
       cerr << "Carrying out dot products for each rotor: " << endl;
       if ( r1.size() != r2.size() ) {
            return 0.0; // mismatch!
       }
       
       mod1 = sqrt( dotdot ( r1 ) );
       mod2 = sqrt( dotdot ( r2 ) );

       double dotter = 0.0;
       double result = 0.0;

     for ( int i = 0 ; i < r1.size() ; i++ ){
         dotter = r1.at(i).dot( r2.at(i) ); // a dot b
         
         cerr << "Poly " << i+1 << " dot: " << dotter << endl;
         
         result += dotter ; // add sum of squares to total
     }       
       
     return result;
}


double dotdot (vector< Vector > &r1 ){
     double result = 0.0;
     double dotter = 0.0;
     
     for ( int i = 0 ; i < r1.size() ; i++ ){
         dotter = r1.at(i).sq() ; // selfdot
         result += dotter ; // add sum of squares to total
     }
     return result;
}

//
//readcom pulls rotors out of the given file
//
bool readcom( string filename, vector< Vector > &rot){
      //outfile << "ID M2-before M2-after bend stretch Rx Ry Rz ModR InDegrees" << endl; //reminder!
     cerr << "Seeking to read rotors from " << filename << endl;
     rot.clear(); // initialise
     Vector vecin; // for reading
     bool pastheader = false;
     int counter = 0; // how many lines read

     string linein;
     string foo; //use as buffer
     vector< string > tokensin;
     
     ifstream comfile( filename.c_str() );
     if ( comfile.fail() ){
        cerr << "Problem! Cannot read file " << filename << endl;
        return false;          
     }
     
     while ( !comfile.eof() ){
            getline(comfile, foo);
            linein = to_lower(foo);// lose case for simplicity
            //cerr << linein << endl;
            tokensin = chop(linein, " \t");
            if (tokensin.size() == 0 ) break; //blank line? must be end of read
            
            counter++ ;// increment line count
            if ( counter == 1 ){
               //that must be the header
               cerr << foo << endl; // just so we can see
            }
            else{
                 //actual line
                 if ( tokensin.size() != 10 ){
                      cerr << "Ending on read of line " << counter << " : " << foo << endl;
                      break; // wrong length of line, end read
                 }
                 double x,y,z ; // components
                 x = atof ( tokensin.at(5).c_str() );
                 y = atof ( tokensin.at(6).c_str());
                 z = atof ( tokensin.at(7).c_str());
                 vecin = Vector( x,y,z );
                 rot.push_back ( vecin ) ;
                 cerr << "Line " << counter << " : vector " << vecin << endl;
            }
            
                             
     }
     if ( counter < 2 ){
        cerr << "No vector data read " << endl;
        return false;
     }
     return true; // successful read       
}


