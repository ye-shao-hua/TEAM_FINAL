#include "../include/fem/hyper_cube.h"

namespace fem{

template<int dim, typename NumberType> HyperCube<dim, NumberType>::HyperCube(
        std::initializer_list<std::initializer_list<NumberType>> list) {
    for(int i = 0; i < list.size(); ++i) {
        _lineXD[i] = {list.begin()[i].begin()[0], list.begin()[i].begin()[1]};
    }
}

template<int dim, typename NumberType> 
const std::pair<NumberType, NumberType>& HyperCube<dim, NumberType>::line(int dimension) const {
    return _lineXD[dimension];
}


template class HyperCube<1,double>;
template class HyperCube<1,float>;
template class HyperCube<2,double>;
template class HyperCube<2,float>;
template class HyperCube<3,double>;
template class HyperCube<3,float>;

}

