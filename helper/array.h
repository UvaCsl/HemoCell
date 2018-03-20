/*
This file is part of the HemoCell library

HemoCell is developed and maintained by the Computational Science Lab 
in the University of Amsterdam. Any questions or remarks regarding this library 
can be sent to: info@hemocell.eu

When using the HemoCell library in scientific work please cite the
corresponding paper: https://doi.org/10.3389/fphys.2017.00563

The HemoCell library is free software: you can redistribute it and/or
modify it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

The library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef HEMO_ARRAY_H
#define HEMO_ARRAY_H

namespace hemo {

  
  template<typename _Tp, std::size_t _Nm>
  struct Array : public std::array<_Tp,_Nm> {
    Array() {}
    Array(const std::initializer_list<_Tp> & il){
      std::copy(il.begin(),il.end(), this->_M_elems);
    }
    Array(const plb::Array<_Tp,_Nm> & ext) {
      for (std::size_t i = 0 ; i < _Nm ; i++) {
        (*this)[i] = ext[i];
      }
    } 
    
    inline Array<_Tp, _Nm> &
    operator+=(const hemo::Array<_Tp, _Nm>& __two) {
      for(std::size_t i = 0 ; i < _Nm ; i++ ) {
        (*this)[i] += __two[i];
      }
      return *this;
    }
    
    inline Array<_Tp, _Nm> &
    operator-=(const hemo::Array<_Tp, _Nm>& __two) {
      for(std::size_t i = 0 ; i < _Nm ; i++ ) {
        (*this)[i] -= __two[i];
      }
      return *this;
    }
    
    inline Array<_Tp, _Nm> &
    operator+=(const plb::Array<_Tp, _Nm>& __two) {
      for(std::size_t i = 0 ; i < _Nm ; i++ ) {
        (*this)[i] += __two[i];
      }
      return *this;
    }
        
    inline Array<_Tp, _Nm> &
    operator*=(const _Tp & mul) {
      for(std::size_t i = 0 ; i < _Nm ; i++ ) {
        (*this)[i] *= mul;
      }
      return *this;
    }
    
    inline Array<_Tp, _Nm> &
    operator/=(const _Tp & div) {
      for(std::size_t i = 0 ; i < _Nm ; i++ ) {
        (*this)[i] /= div;
      }
      return *this;
    }
    
    
    inline Array<_Tp, _Nm> &
    operator=(const std::initializer_list<_Tp> & c) {
      std::copy(c.begin(),c.end(), this->_M_elems);
      return *this;
    } 
    inline void resetToZero() {
      this->fill(0);
    }

  };
  
  template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator+(const Array<_Tp, _Nm> & one, const Array<_Tp, _Nm> & two) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]+two[i];
    }
    return ret;
  }
  
    template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator+(const Array<_Tp, _Nm> & one, const plb::Array<_Tp, _Nm> & two) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]+two[i];
    }
    return ret;
  }
  template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator+(const Array<_Tp, _Nm> & one, const _Tp & two) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]+two;
    }
    return ret;
  }
  template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator-(const Array<_Tp, _Nm> & one, const Array<_Tp, _Nm> & two) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]-two[i];
    }
    return ret;
  }
  template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator-(const Array<_Tp, _Nm> & one) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = -one[i];
    }
    return ret;
  }

  template<typename _Tp, std::size_t _Nm, typename _Tp2>
  inline Array<_Tp, _Nm> operator/(const Array<_Tp, _Nm> & one, const _Tp2 div) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]/div;
    }
    return ret;
  }
  
  template<typename _Tp, std::size_t _Nm, typename _Tp2>
  inline Array<_Tp, _Nm> operator*(const Array<_Tp, _Nm> & one, const _Tp2 mul) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]*mul;
    }
    return ret;
  }
  
  template<typename _Tp, std::size_t _Nm, typename _Tp2>
  inline Array<_Tp, _Nm> operator*(const _Tp2 mul, const Array<_Tp, _Nm> & one) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]*mul;
    }
    return ret;
  }
  
  template<typename _Tp>
  inline void crossProduct (const Array<_Tp, 3>& one, const Array<_Tp,3> & two, Array<_Tp,3> & result) {
    result[0] = one[1]*two[2] - one[2]*two[1];
    result[1] = one[2]*two[0] - one[0]*two[2];
    result[2] = one[0]*two[1] - one[1]*two[0];
  };
   
  template<typename _Tp>
  inline Array<_Tp,3> crossProduct (const Array<_Tp, 3>& one, const Array<_Tp,3> & two) {
    Array<_Tp,3> result;
    crossProduct(one, two, result);
    return result;
  };
  
  template<typename _Tp, std::size_t _Nm>
  inline _Tp dot(const Array<_Tp, _Nm> & one, const Array<_Tp, _Nm> & two) {
    _Tp ret = 0;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret += one[i]*two[i];
    }
    return ret;
  };
  
  template<typename _Tp, std::size_t _Nm>
  inline _Tp norm(const Array<_Tp, _Nm> & one) {
    _Tp ret = 0;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret += (one[i])*(one[i]);
    }
    return std::sqrt(ret);
  };
  
  template<typename _Tp>
  _Tp angleBetweenVectors(const Array<_Tp,3> & one, const Array<_Tp,3> & two)
  {
    Array<_Tp,3> cross;
    crossProduct(one, two, cross); 
    return std::atan2(norm(cross), dot(one,two));
  }
 
  template<typename _Tp>
  _Tp computeTriangleArea(const Array<_Tp,3> & v0, const Array<_Tp,3> & v1, const Array<_Tp,3> & v2)
  {
    Array<_Tp,3> e1 = v1 - v0;
    Array<_Tp,3> e2 = v2 - v0;
    Array<_Tp,3> cross;
    crossProduct(e1, e2, cross);
    
    return (_Tp) 0.5 * norm(cross);
  }
  
  template<typename _Tp>
  void computeTriangleAreaAndUnitNormal(const Array<_Tp,3> & v0, const Array<_Tp,3> & v1, const Array<_Tp,3> & v2, _Tp & area, Array<_Tp,3>& unitNormal)
  {
      Array<_Tp,3> e01 = v1 - v0;
      Array<_Tp,3> e02 = v2 - v0;

      crossProduct(e01, e02, unitNormal);
      _Tp normN = norm(unitNormal);
      if (normN != 0.0) {
          area = (_Tp) 0.5 * normN;
          unitNormal /= normN;
      } else {
          area = (_Tp) 0;
          unitNormal.resetToZero();
      }
  }
  
  template<typename _Tp>
  Array<_Tp,3> computeTriangleNormal(const Array<_Tp,3> & v0, const Array<_Tp,3> & v1, const Array<_Tp,3> & v2, bool isAreaWeighted)
  {
      Array<_Tp,3> e01 = v1 - v0;
      Array<_Tp,3> e02 = v2 - v0;

      Array<_Tp,3> n;
      crossProduct(e01, e02, n);
      if (!isAreaWeighted) {
          _Tp normN = norm(n);
          if (normN != 0) {
              n /= normN;
          } else {
              n.resetToZero();
          }
      }
      return n;
  }
  
  template<typename _Tp>
  _Tp computeCotangentFromVectors (const Array<_Tp,3> & a, const Array<_Tp,3> & b) {
    const _Tp adj_norm = dot(a,b)/norm(b);
    //Pythagoras to get opp side length
    const _Tp hyp_norm = norm(a);
    const _Tp opp_norm = std::sqrt(std::pow(hyp_norm,2)-std::pow(adj_norm,2));
    return adj_norm/opp_norm;
  }
  
  template<typename _Tp>
  void computeLengthsPythagoras (const Array<_Tp,3> & a, const Array<_Tp,3> & b, _Tp & a_l, _Tp & b_l, _Tp & c_l) {
    a_l = dot(a,b)/norm(b);
    //Pythagoras to get opp side length
    c_l = norm(a);
    b_l = std::sqrt(std::pow(c_l,2)-std::pow(a_l,2));
  }  
}

#endif
