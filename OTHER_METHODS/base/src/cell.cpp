#include "./../include/fem/cell.h"
#include "./../include/fem/base_function.h"
#include "./../include/fem/point.h"
#include "./../include/fem/hyper_cube.h"
namespace fem{
    
    template<int dim,typename NumberType>
        Cell<dim,NumberType>::Cell(const HyperCube<dim,NumberType>& size):_size(size){}

    template<int dim,typename NumberType>
    NumberType Cell<dim,NumberType>::operator()(int id,const Point<dim,NumberType> &position){
        return _base->operator()(id,position);
    }

    template<int dim,typename NumberType>
    const HyperCube<dim,NumberType> &Cell<dim,NumberType>::get_size(){
        return _size;
    }

    

template class Cell<1,double>; 
template class Cell<1,float>; 
template class Cell<2,double>; 
template class Cell<2,float>; 
template class Cell<3,double>; 
template class Cell<3,float>; 
}
