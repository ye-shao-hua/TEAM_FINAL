#include "./../include/fem/fe_q.h"
#include "./../include/fem/cell.h"
#include "./../include/fem/hyper_cube.h"
namespace fem {
    template<int dim,typename NumberType>
        FE_Q<dim,NumberType>::FE_Q(int max_degree):_max_degree(max_degree){}

    template<int dim,typename NumberType>
        void FE_Q<dim,NumberType>::init(int length){
            _base_functions.reserve(length);
            for(auto i=0;i<length; ++i){
                _base_functions.emplace_back(_max_degree);
            }
        }
    
    template<int dim,typename NumberType>
        void FE_Q<dim,NumberType>::distribute(Cell<dim,NumberType> &cell){
            static int cnt=0;
            cell._base=_base_functions.begin()+cnt++;
        }

    template class FE_Q<1,double>;
    template class FE_Q<1,float>;
    template class FE_Q<2,double>;
    template class FE_Q<2,float>;
    template class FE_Q<3,double>;
    template class FE_Q<3,float>;
}
