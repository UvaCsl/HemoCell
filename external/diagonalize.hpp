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

// header of 3-vectorDiag

class vectorDiag_exception : public exception
{ public:
    vectorDiag_exception( const string & );
    ~vectorDiag_exception() throw();
  const char * what() const throw();
  protected:
    string error_description;
};

class matrixDiag;

class vectorDiag
{ double the_data[ 3 ];
 public:
  friend class matrixDiag;
  vectorDiag();
  vectorDiag( const vectorDiag & );
  ~vectorDiag();
  double & operator[]( unsigned int );
  const double & operator[]( unsigned int ) const;
  vectorDiag & set_zero();
  vectorDiag & operator=( const vectorDiag & );
  vectorDiag & operator+=( const vectorDiag & );
  vectorDiag & operator-=( const vectorDiag & );
};

istream & operator>>( istream &, vectorDiag & );
ostream & operator<<( ostream &, const vectorDiag & );

// header of 3x3 matrixDiag

class matrixDiag_exception : public exception
{ public:
    matrixDiag_exception( const string & );
    ~matrixDiag_exception() throw();
  const char * what() const throw();
  protected:
    string error_description;
};

class matrixDiag
{ double the_data[ 9 ];
 public:
  friend class vectorDiag;
  matrixDiag();
  matrixDiag( const matrixDiag & );
  ~matrixDiag();
  class array_ref
  { friend class matrixDiag;
    matrixDiag * the_ref;
    unsigned int the_pos;
    array_ref( matrixDiag *, unsigned int );
   public:
    ~array_ref();
    double & operator[]( unsigned int );
    const double & operator[]( unsigned int ) const;
  };
  friend class array_ref;
  array_ref operator[]( unsigned int );
  const array_ref operator[]( unsigned int ) const;
  double norm() const;
  matrixDiag & set_zero();
  matrixDiag & transpose();
  matrixDiag & set_rotation( unsigned int, double );
  matrixDiag & set_diagonal( const vectorDiag & );
  matrixDiag & sym_eigen( vectorDiag &, vectorDiag & );
  matrixDiag & operator=( const matrixDiag & );
  matrixDiag & operator+=( const matrixDiag & );
  matrixDiag & operator-=( const matrixDiag & );
  matrixDiag & operator*=( const matrixDiag & );
};

istream & operator>>( istream &, matrixDiag & );
ostream & operator<<( ostream &, const matrixDiag & );

#endif
