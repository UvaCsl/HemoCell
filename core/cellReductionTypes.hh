#ifndef CELL_REDUCTION_TYPES_HH
#define CELL_REDUCTION_TYPES_HH

#include "cellReductionTypes.h"


/*******************************************************************
 *                      1D Subscribe Reductions                    *
 *******************************************************************/
template<typename T>
plint BlockStatisticsCCR<T>::subscribeReduction1D(plint reductionType) {
    if (0 == reductionType)      { return subscribeSum(); }
    else if (1 == reductionType) { return subscribeAverage(); }
    else if (2 == reductionType) { return subscribeMax(); }
    else if (3 == reductionType) { return subscribeMax(); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { return subscribeAverage(); } // Std is essentially an average
    else { return -1; }
}

/*******************************************************************
 *                         1D Gather Reductions                    *
 *******************************************************************/

template<typename T>
void BlockStatisticsCCR<T>::gatherReduction1D(plint reductionType, plint qBin, T value) {
    if (0 == reductionType)      { gatherSum(qBin, value); }
    else if (1 == reductionType) { gatherAverage(qBin, value); }
    else if (2 == reductionType) { gatherMax(qBin, value); }
    else if (3 == reductionType) { gatherMax(qBin, -value); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { gatherAverage(qBin, value); } // Std is essentially an average
}


/*******************************************************************
 *                          1D Get Reductions                      *
 *******************************************************************/

template<typename T>
T BlockStatisticsCCR<T>::getReduction1D(plint reductionType, plint qBin) {
    if (0 == reductionType)      { return getSum(qBin); }
    else if (1 == reductionType) { return getAverage(qBin); }
    else if (2 == reductionType) { return getMax(qBin); }
    else if (3 == reductionType) { return -getMax(qBin); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { return getAverage(qBin); } // Std is essentially an average
    return -1;
}



/*
 * 3D and ND Reduction operations
 */
// subscribeReduction
template<typename T>
Array<plint,3> BlockStatisticsCCR<T>::subscribeReduction3D(plint reductionType) {
    plint x = subscribeReduction1D(reductionType);
    plint y = subscribeReduction1D(reductionType);
    plint z = subscribeReduction1D(reductionType);
    return Array<plint,3>(x, y, z);
}

template<typename T>
std::vector<plint> BlockStatisticsCCR<T>::subscribeReductionND(plint reductionType, plint dimensions) {
    std::vector<plint> ret;
    for (plint i = 0; i < dimensions; ++i) { ret.push_back( subscribeReduction1D(reductionType) ); }
    return ret;
}

// gatherReduction
template<typename T>
void BlockStatisticsCCR<T>::gatherReduction3D(plint reductionType, Array<plint,3> qBin, Array<T,3> value) {
    for (int i = 0; i < 3; ++i) { gatherReduction1D(reductionType, qBin[i], value[i]); }
}

template<typename T>
void BlockStatisticsCCR<T>::gatherReductionND(plint reductionType, std::vector<plint> qBin, std::vector<T> value) {
    for (pluint i = 0; i < qBin.size(); ++i) { gatherReduction1D(reductionType, qBin[i], value[i]); }
}

// getReduction
template<typename T>
Array<T,3> BlockStatisticsCCR<T>::getReduction3D(plint reductionType, Array<plint,3> qBin) {
    T x = getReduction1D(reductionType, qBin[0]);
    T y = getReduction1D(reductionType, qBin[1]);
    T z = getReduction1D(reductionType, qBin[2]);
    return Array<T,3>(x, y, z);
}

template<typename T>
std::vector<T> BlockStatisticsCCR<T>::getReductionND(plint reductionType, std::vector<plint> qBin) {
    std::vector<T> ret;
    plint dimensions = qBin.size();
    for (plint i = 0; i < dimensions; ++i) { ret.push_back( getReduction1D(reductionType, qBin[i]) ); }
    return ret;
}



#endif
