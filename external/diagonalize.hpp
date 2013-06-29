/*
   DISCLAIMER OF WARRANTY
   The author makes no warranties, express or implied, that this program
   is free of error, or consistent with any particular standard of
   merchantability, or that it will meet your requirements for any
   particular application. It should not be relied on for solving a problem
   whose incorrect solution could result in injury to a person or loss
   of property. If you do use the program in such a manner, it is at your own
   risk. The author disclaims all liability for direct or consequential
   damages resulting from your use of the program.
*/
//---------------------------------------------------------------------------
#ifndef DiagonalizeHpp
#define DiagonalizeHpp
//---------------------------------------------------------------------------

#include <stdexcept>
#include <string>

using namespace std;

// header of 3-vector

class vector_exception : public exception
{ public:
    vector_exception( const string & );
    ~vector_exception() throw();
  const char * what() const throw();
  protected:
    string error_description;
};

class matrix;

class vector
{ double the_data[ 3 ];
 public:
  friend class matrix;
  vector();
  vector( const vector & );
  ~vector();
  double & operator[]( unsigned int );
  const double & operator[]( unsigned int ) const;
  vector & set_zero();
  vector & operator=( const vector & );
  vector & operator+=( const vector & );
  vector & operator-=( const vector & );
};

istream & operator>>( istream &, vector & );
ostream & operator<<( ostream &, const vector & );

// header of 3x3 matrix

class matrix_exception : public exception
{ public:
    matrix_exception( const string & );
    ~matrix_exception() throw();
  const char * what() const throw();
  protected:
    string error_description;
};

class matrix
{ double the_data[ 9 ];
 public:
  friend class vector;
  matrix();
  matrix( const matrix & );
  ~matrix();
  class array_ref
  { friend class matrix;
    matrix * the_ref;
    unsigned int the_pos;
    array_ref( matrix *, unsigned int );
   public:
    ~array_ref();
    double & operator[]( unsigned int );
    const double & operator[]( unsigned int ) const;
  };
  friend class array_ref;
  array_ref operator[]( unsigned int );
  const array_ref operator[]( unsigned int ) const;
  double norm() const;
  matrix & set_zero();
  matrix & transpose();
  matrix & set_rotation( unsigned int, double );
  matrix & set_diagonal( const vector & );
  matrix & sym_eigen( vector &, vector & );
  matrix & operator=( const matrix & );
  matrix & operator+=( const matrix & );
  matrix & operator-=( const matrix & );
  matrix & operator*=( const matrix & );
};

istream & operator>>( istream &, matrix & );
ostream & operator<<( ostream &, const matrix & );

#endif
