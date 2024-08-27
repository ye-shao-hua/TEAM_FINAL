#include "../include/fem/point.h"

namespace fem{
template<int dim,typename NumberType> 
    Point<dim,NumberType>::Point(const std::array<NumberType,dim> &data):_data(data){}

template<int dim,typename NumberType>
    NumberType& Point<dim,NumberType>::
    get(int index) {
        return _data[index];
    }

template<int dim,typename NumberType>
    const NumberType& Point<dim,NumberType>::
    get(int index) const {
        return _data[index];
    }

template class Point<1,double>;
template class Point<2,double>;
template class Point<3,double>;
template class Point<1,float>;
template class Point<2,float>;
template class Point<3,float>;

}

