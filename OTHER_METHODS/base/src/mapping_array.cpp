#include "./../include/fem/mapping_array.h"
#include "../include/fem/point.h"
#include "../include/fem/full_matrix.h"
namespace fem{
    template<int dim,typename NumberType>
        MappingArray<dim,NumberType>::MappingArray(Point<dim,NumberType> *point){
            //按顺序传入四个点
            //先判断是矩形还是梯形
            if((point[0].get(0)==point[3].get(0)) && (point[1].get(0)==point[2].get(0))){
                //矩形
            }else{
                //梯形
            }
        }
}
