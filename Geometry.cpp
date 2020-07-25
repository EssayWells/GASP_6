#include "Geometry.h"

//
//printbond returns a bondline entry in printable form
//
string Bondline::printbond ( ){
       stringstream result;
       result << id1+1 << " " << id2+1 << "   " << bars << "  " << ele1 << " " << ele2;
       return result.str();
       //the +1 is because internally atoms are numbered from 0 upwards
}



//
//readbonds pulls bond data out of the given file and into the bondline array
//
bool readbonds( string filename, vector< Bondline > &bl ){
     cerr << "Seeking to read bonding from " << filename << endl;
     bl.clear();

     string linein;
     string foo; //use as buffer
     vector< string > tokensin;
     
     ifstream bondfile( filename.c_str() );
     if ( bondfile.fail() ) {
        cerr << "Problem! Cannot read file " << filename << endl;
        return false;     
     }

     while( getline(bondfile, foo) ){
            linein = to_lower(foo);// lose case for simplicity
            //cerr << linein << endl;
            tokensin = chop(linein, " \t");
            if (tokensin.size() == 0 ) continue; //blank line?
            
            if ( tokensin.size() < 5 ) {
                 cerr << "Error: malformed line " << foo << " in bond file " << endl;
                 return false;
            }
            
            Bondline dummy;
            dummy.id1 = atoi( tokensin.at(0).c_str() ) -1;//offset!
            dummy.id2 = atoi( tokensin.at(1).c_str() ) -1;//c++ 0 up, humans 1 up
            dummy.bars = atoi( tokensin.at(2).c_str() );
            dummy.ele1 = tokensin.at(3);
            dummy.ele2 = tokensin.at(4);
            
            bl.push_back( dummy );
     }
     bondfile.close();
     return true; // all done!
}


//
//writebonds exports the bonding of the structure into a file
//note "5 bars" = rotatable covalent bond, conventional;
//"6 bars" = nonrotatable covalent bond, connects sp2 carbons or amides
//printbond will +1 the atom ids for output
void writebonding( string filename, vector< Bondline > &bl ){
     int nb = bl.size();
     cerr << "Writing bonding to " << filename << endl;
     ofstream bondfile( filename.c_str() );     
     for ( int i = 0 ; i < nb ; i++ ){
         bondfile << bl.at(i).printbond() << endl;
     }
     
     bondfile.close();
     return;
}


