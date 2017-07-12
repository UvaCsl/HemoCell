/* 
 * File:   palabos_array_interface.h
 * 
 * This file contains the (minimal) interface to plb::Array, only used to make
 * interaction with Palabos possible (without even linking etc). This to
 * minimize GPL effects
 * 
 */



#ifndef PALABOS_ARRAY_INTERFACE_H
#define PALABOS_ARRAY_INTERFACE_H

namespace plb {

/** 
 * A interface to the PLB Array class.
 */
template<typename T, size_t size>
class Array {
public:
    T& operator[](size_t index) {
        return data[index];
    }
    T const& operator[](size_t index) const {
        return data[index];
    }
private:
    T data[size];
};

}

#endif /* PALABOS_ARRAY_INTERFACE_H */

