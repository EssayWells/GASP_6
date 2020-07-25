#include "basic_utils.h"

//this utils package includes mt19337 with associated disclaimers, see below

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          
   Copyright (C) 2005, Mutsuo Saito,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdio.h>
//#include "mt19937ar.h"

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */




////////////////////////////////////////////////////////////////////////////////
// Description:
//   Convert a string to all lower-case letters.
// Parameters:
//   str: Input string.
////////////////////////////////////////////////////////////////////////////////
string to_lower( string &str ){

  for( unsigned int a = 0; a < str.length(); a++ ){
    if( str[a] >= 'A' && str[a] <= 'Z' )
      str[a] += 32;
  }
  
  return(str);
}

////////////////////////////////////////////////////////////////////////////////
// Description:
//   Convert a string to all upper-case letters.
// Parameters:
//   str: Input string.
////////////////////////////////////////////////////////////////////////////////
string to_upper( string &str ){

  for( unsigned int a = 0; a < str.length(); a++ ){
    if( str[a] >= 'a' && str[a] <= 'z' )
      str[a] -= 32;
  }
  
  return(str);
}

////////////////////////////////////////////////////////////////////////////////
// Description: Generic string tokenizer. The function takes as arguments the
//   string to be parsed and an optional string containing the delimiters. By 
//   default, the function will assume white-space characters as delimiters. The
//   function returns a vector of strings containing the tokens parsed from the
//   input string. NOTE: If your input string contains newline characters, these
//   will be treated as a delimiter. The string will continue past intervening 
//   newlines until it reaches the end of the string, as determined by string::npos.
////////////////////////////////////////////////////////////////////////////////
vector<string> chop( string &input, string delimiters ){

  // Add the newline character to the list of tokens, just in case it 
  // wasn't included.
  //////////////////////////////////////////////////////////////////////
  if( delimiters.find("\n") == string::npos )
    delimiters += "\n";

  vector<string> tokens;

  size_t start = input.find_first_not_of(delimiters,0);
  size_t end   = input.find_first_of(delimiters,start);

  if( start == string::npos )
    return( tokens );

  while( end != string::npos ){
    tokens.push_back( input.substr(start, end-start) );
    start = input.find_first_not_of( delimiters, end );
    end   = input.find_first_of( delimiters, start );
  }

  if( start != string::npos )
    tokens.push_back( input.substr(start, end-start) );

  return( tokens );
}


////////////////////////////////////////////////////////////////////////////////
// Description: Remove leading and trailing white-space characters from a 
//   string. 
////////////////////////////////////////////////////////////////////////////////
void chomp( string &str ){

  if ( str == "" ) return; //empty strings

  int start = str.find_first_not_of(" \t");
  int end   = str.find_last_not_of(" \t\n");

  //cerr << str << " [" << str.substr(start, end-start+1) << "] " << start << " " << end << endl;
  str = str.substr( start, end-start+1 );
}

void initialiseRandom(unsigned long s){
	init_genrand(s);
}
double getRandom(){
	double r;
	r = genrand_real1();
	return r;
}

bool isin ( int id, vector< int > &list ){
     for ( int i = 0 ; i < list.size(); i++ ){
         if ( list.at(i) == id ) return true;
     }
     return false;
}

string lowerCaseVersion( string str ){
	string result;
	result = to_lower( str ); //passed value not reference, shouldn't alter the argument
	return result;
}

bool almostEqual( double i, double j, double margin /* = 1.e-4 */ ){
	//returns true if two doubles differ by only a small amount
	//found in testing:
	//two doubles, both equal to "90", differ by more than 1e-6!
	//that's terrible.
	
	//double vsmall = 1.e-4; //is now a default value
	//cerr << "Values: " << i << " , " << j << endl;
	
	double d = i - j;
	if ( d < 0. ) d *= -1.; //flip!
	//cerr << "Difference: " << d << endl;
	if ( d < margin ) return true;
	//cerr << "Beep beep." << endl;
	return false;
}



