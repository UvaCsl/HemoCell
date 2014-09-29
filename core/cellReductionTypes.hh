#ifndef CELL_REDUCTION_TYPES_HH
#define CELL_REDUCTION_TYPES_HH

#include "cellReductionTypes.h"



template<typename T>
void BlockStatisticsCCR<T>::subscribe(plint ccrId) {
    plint reductionType = getReductionType(ccrId);
    plint dim = getReductionDimension(ccrId);
    pluint qBin;
    for (int d = 0; d < dim; ++d) {
        if (0 == reductionType)      { qBin = subscribeSum(); }
        else if (1 == reductionType) { qBin = subscribeAverage(); }
        else if (2 == reductionType) { qBin = subscribeMax(); }
        else if (3 == reductionType) { qBin = subscribeMax(); } // Min can be calculated from the inverse Max
        else if (4 == reductionType) { qBin = subscribeAverage(); } // Std is essentially an average
        else { qBin = -1; }
    }
    ccrIdToBin[ccrId] = qBin - (dim-1);
}

template<typename T>
void BlockStatisticsCCR<T>::gather(plint ccrId, T value, pluint qBin) {
    plint reductionType = getReductionType(ccrId);
    if (0 == reductionType)      { gatherSum(qBin, value); }
    else if (1 == reductionType) { gatherAverage(qBin, value); }
    else if (2 == reductionType) { gatherMax(qBin, value); }
    else if (3 == reductionType) { gatherMax(qBin, -value); } // Min can be calculated from the inverse Max
    else if (4 == reductionType) { gatherAverage(qBin, value); } // Std is essentially an average
}


template<typename T>
void BlockStatisticsCCR<T>::gather(plint ccrId, T value) {
    pluint qBin = ccrIdToBin[ccrId];
    gather(ccrId, value, qBin);
}


template<typename T>
void BlockStatisticsCCR<T>::gather(plint ccrId, Array<T,3> const& value) {
    plint dim = getReductionDimension(ccrId);
    pluint qBin;
    for (int d = 0; d < dim; ++d) {
        qBin = ccrIdToBin[ccrId] + d;
        gather(ccrId, value[d], qBin);
    }
}


template<typename T>
void BlockStatisticsCCR<T>::gather(plint ccrId, std::vector<T> const& value) {
    plint dim = getReductionDimension(ccrId);
    pluint qBin;
    for (int d = 0; d < dim; ++d) {
        qBin = ccrIdToBin[ccrId] + d;
        gather(ccrId, value[d], qBin);
    }
}

template<typename T>
void BlockStatisticsCCR<T>::get(plint ccrId, T & value, pluint qBin) {
    gather(ccrId, value, qBin);
}

template<typename T>
void BlockStatisticsCCR<T>::get(plint ccrId, T & value) {
    pluint qBin = ccrIdToBin[ccrId];
    get(ccrId, value, qBin);
}


template<typename T>
void BlockStatisticsCCR<T>::get(plint ccrId, Array<T,3> & value) {
    plint dim = getReductionDimension(ccrId);
    pluint qBin;
    for (int d = 0; d < dim; ++d) {
        qBin = ccrIdToBin[ccrId] + d;
        get(ccrId, value[d], qBin);
    }
}


template<typename T>
void BlockStatisticsCCR<T>::get(plint ccrId, std::vector<T> & value) {
    plint dim = getReductionDimension(ccrId);
    pluint qBin;
    value.resize(dim);
    for (int d = 0; d < dim; ++d) {
        qBin = ccrIdToBin[ccrId] + d;
        gather(ccrId, value[d], qBin);
    }
}





#endif
