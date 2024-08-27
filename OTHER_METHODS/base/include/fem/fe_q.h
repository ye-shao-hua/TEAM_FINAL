#pragma once
#include "./base_function.h"
#include <vector>
namespace fem{

    template<int dim,typename NumberType> class Cell;

    template<int dim,typename NumberType>
    class FE_Q {
        public:
        FE_Q(int max_degree);
        void init(int length);
        void distribute(Cell<dim,NumberType> &cell);
        private:
        std::vector<BaseFunctionGroupXD<dim,NumberType>> _base_functions;
        int _max_degree;
    };
}
