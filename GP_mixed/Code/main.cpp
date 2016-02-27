#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iterator>
#include <armadillo>
#include <cmath> //for erfc function and math constants
#include "eOfZ.h"
#include "eOfZTZ.h"
#include "logH.h"
#include "lowerBound.h"
#include "sparseGPVBProbit.h"
#include "gaussianCDF.h"
#include "mvnRandomGenerate.h"
using namespace std;


template <class T>
class csv_istream_iterator: public iterator<input_iterator_tag, T>
{
    istream * _input;
    char _delim;
    string _value;
public:
    csv_istream_iterator( char delim = ',' ): _input( 0 ), _delim( delim ) {}
    csv_istream_iterator( istream & in, char delim = ',' ): _input( &in ), _delim( delim ) { ++*this; }

    const T operator *() const {
        istringstream ss( _value ); 
        T value;
        ss >> value;
        return value;
    }

    istream & operator ++() {
        if( !( getline( *_input, _value, _delim ) ) )
        {
            _input = 0;
        }
        return *_input;
    }

    bool operator !=( const csv_istream_iterator & rhs ) const {
        return _input != rhs._input;
    }
};

int main() {
  // for simulation
  // The true regression model is P(y = 1|X, β) = Φ(0.3 + (0.216 * x1) + (0.1264 * x2) + (2.375 * x3) + a lot of random effects)
  vector<float> linear;

  ifstream fin( "matrix.csv" );
  if( fin )
  {
      copy( csv_istream_iterator<float>( fin ),
            csv_istream_iterator<float>(),
            insert_iterator< vector<float> >( linear, linear.begin() ) );

      fin.close();
  }

  // linear now contains all floating point values in the comma-delimited
  // data file.  now dump them to the screen to verify:
  // copy( linear.begin(), linear.end(), ostream_iterator<float>( cout, " " ) );

  arma::mat tempMat = arma::conv_to<arma::mat>::from(linear);
  int n = tempMat.n_rows;
  arma::mat temp(n-1, 1);
  for (int i = 1; i < n; i++) {
    temp(i-1, 0) = tempMat(i, 0);
  }
  temp.reshape(9, 315);
  arma::mat X(315, 9);
  X = temp.t();
  

  arma::vec beta = arma::randn<arma::vec>(9);
  arma::vec constructY = X * beta;
  arma::colvec y(315);
  for (int eachRow = 0; eachRow < 315; eachRow++) {
    if (constructY(eachRow) >= 0.0) {
      y(eachRow) = 1.0;
    } else {
      y(eachRow) = 0.0;
    }
  }

  


  return 0;
}