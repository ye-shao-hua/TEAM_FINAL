#include "../include/fem/mapping_q.h"
#include "../include/fem/full_matrix.h"
#include <iostream>
#include "../include/fem/point.h"
#include <cmath>
#include <array>

namespace fem {

template <int dim, typename NumberType>
MappingQ<dim, NumberType>::MappingQ(int degree) {
  for (int i = 0; i < degree; ++i)
    _data.emplace_back(dim, i + 1);
}

template <int dim, typename NumberType>
void MappingQ<dim, NumberType>::set_jaccobi(int id, const FullMatrix &matrix) {
  _data[id] = matrix;
}


template <int dim, typename NumberType>
Point<dim,NumberType> MappingQ<dim,NumberType>::mapping(Point<dim,NumberType> point){
    //NumberType first = pair1.first;
    //NumberType second = pair1.second;
    Point<dim,NumberType> result{std::array<NumberType,dim>{}};
    if(dim==2){
        for(auto i:_data){
            int cols=i.cols();
            FullMatrix full1(cols,1);
            for(auto j=0;j<cols;++j){
                full1(j,0)=std::pow(point.get(1),j)*std::pow(point.get(0),cols-1-j);
            }
            result.get(0)+=(i*full1)(0,0);
            result.get(1)+=(i*full1)(1,0);
        
        }
    }else if(dim!=1){
        std::cout<<"error for mapping\n";
    }else{
        for(auto i=0;i<_data.size();++i){
            int cols=_data[i].cols();
            FullMatrix full1(cols,1);
            full1(0,0)=std::pow(point.get(0),i);
            result.get(0)+=(_data[i]*full1)(0,0);
        }

    }
    return result;
}

template class MappingQ<1, double>;
template class MappingQ<1, float>;
template class MappingQ<2, double>;
template class MappingQ<2, float>;
template class MappingQ<3, double>;
template class MappingQ<3, float>;

} // namespace fem
