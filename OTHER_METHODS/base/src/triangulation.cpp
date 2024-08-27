#include "./../include/fem/triangulation.h"
#include "./../include/fem/hyper_cube.h"
#include "./../include/fem/cell.h"
#include <iostream>
namespace fem{
    template<int dim, typename NumberType>
        Triangulation<dim,NumberType>::Triangulation(HyperCube<dim,NumberType> &hy){
            _data.push_back(Cell<dim,NumberType>{hy});
        }
        
    template class Triangulation<1,double>;
    template class Triangulation<1,float>;
    template class Triangulation<2,double>;
    template class Triangulation<2,float>;
    template class Triangulation<3,double>;
    template class Triangulation<3,float>;
}

