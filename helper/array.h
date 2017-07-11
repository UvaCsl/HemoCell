#ifndef HEMO_ARRAY_H
#define HEMO_ARRAY_H

namespace hemo {

  
  template<typename _Tp, std::size_t _Nm>
  struct Array : public std::array<_Tp,_Nm> {
    Array() {}
    Array(const std::initializer_list<_Tp> & il){
      std::copy(il.begin(),il.end(), this->_M_elems);
    }

    inline Array<_Tp, _Nm> &
    operator+=(const hemo::Array<_Tp, _Nm>& __two) {
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
    operator=(const std::initializer_list<_Tp> & c) {
      std::copy(c.begin(),c.end(), this->_M_elems);
      return *this;
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
  inline Array<_Tp, _Nm> operator-(const Array<_Tp, _Nm> & one, const Array<_Tp, _Nm> & two) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]-two[i];
    }
    return ret;
  }
  
  template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator/(const Array<_Tp, _Nm> & one, const _Tp div) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]/div;
    }
    return ret;
  }
  
  template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator*(const Array<_Tp, _Nm> & one, const _Tp mul) {
    Array<_Tp, _Nm> ret;
    for (std::size_t i = 0 ; i < _Nm ; i++ ) {
      ret[i] = one[i]*mul;
    }
    return ret;
  }
  
  template<typename _Tp, std::size_t _Nm>
  inline Array<_Tp, _Nm> operator*(const _Tp mul, const Array<_Tp, _Nm> & one) {
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
}
#endif
