#include "../include/fem/base_function.h"
#include <cmath>
#include "../include/fem/point.h"

namespace fem {

template <typename NumberType> NumberType AttenuationFunction1D(double x) {
  return std::max(std::min(-x, x) + 1., 0.);
}

template double AttenuationFunction1D<double>(double x);
template float AttenuationFunction1D<float>(double x);

template <typename NumberType>
NumberType BaseFunction1D<NumberType>::operator()(double x) {
  return std::pow(x, _degree);
}
template <typename NumberType>
BaseFunction1D<NumberType>::BaseFunction1D(int degree):_degree(degree){}

template<typename NumberType>
NumberType BaseFunction1D<NumberType>::grad(double x){
    return _degree==0?0:_degree* std::pow(x,_degree-1);
}

template class BaseFunction1D<double>;
template class BaseFunction1D<float>;

template<int dim, typename NumberType> 
BaseFunctionGroupXD<dim, NumberType>::
    BaseFunctionGroupXD(int max_degree){
    for(std::vector<BaseFunction1D<NumberType>> &i:_base_functions){
        for(auto j=0;j<=max_degree;++j){
            i.template emplace_back<BaseFunction1D<NumberType>>(j);
        }
    }    
}

template<int dim, typename NumberType> 
NumberType BaseFunctionGroupXD<dim, NumberType>::
    operator()(int id, const Point<dim, NumberType> &position) {
        std::array<int, dim> indexs;
        for(int i = 0; i < dim; ++i) {
            indexs[i] = id % _base_functions[0].size();
            id /= _base_functions[0].size();
        }
        NumberType value = 1;
        for(int i = 0; i < dim; ++i)
            value *= _base_functions[i][indexs[i]](position.get(i));
        return value;
    }

template<int dim, typename NumberType>
Point<dim, NumberType> BaseFunctionGroupXD<dim, NumberType>::
    grad(int id,const Point<dim,NumberType> &position){
        std::array<int, dim> indexs;
        for(int i = 0; i < dim; ++i) {
            indexs[i] = id % _base_functions[0].size();
            id /= _base_functions[0].size();
        }
        Point<dim,NumberType> result{std::array<NumberType, dim>{}};
        for(int i = 0; i < dim; ++i) {
            result.get(i) = _base_functions[i][indexs[i]].grad(position.get(i));
        }
        return result;
    }

template class BaseFunctionGroupXD<1, double>;
template class BaseFunctionGroupXD<2, double>;
template class BaseFunctionGroupXD<3, double>;
template class BaseFunctionGroupXD<1, float>;
template class BaseFunctionGroupXD<2, float>;
template class BaseFunctionGroupXD<3, float>;

} // namespace fem
