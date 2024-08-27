#pragma once
#include "./point.h"
#include "./base_function.h"
#include <vector>
#include "./mapping_q.h"
namespace fem {
    template <int dim, typename NumberType>
    class HyperCube;

    template<int dim,typename NumberType>
    class Cell{
        template<int dim1,typename NumberType1> friend class FE_Q;
        public:
            Cell(const HyperCube<dim, NumberType>& size);
            NumberType operator()(int id,const Point<dim,NumberType> &position);
            const HyperCube<dim,NumberType> &get_size();
        private:
            HyperCube<dim, NumberType> _size; // 非单位
            //MappingQ<dim,NumberType> _map;
            // MappingQ  
            // _begin_point
            // (x, y) - begin = (diff x, diff y)
            // diff x * dy/dx = 
            // diff y * dy/dy = 
            // (x^, y^)
            // shape_value(id, x^, y^)
            typename std::vector<BaseFunctionGroupXD<dim,NumberType>>::iterator _base;
    };
} 
