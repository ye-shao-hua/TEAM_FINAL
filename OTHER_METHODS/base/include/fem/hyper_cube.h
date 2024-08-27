#pragma once
#include <array>
#include <initializer_list>

namespace fem{

    template<int dim,typename NumberType>
    class HyperCube{
        public:
            HyperCube(std::initializer_list<std::initializer_list<NumberType>> list);
            const std::pair<NumberType, NumberType>& line(int dimension) const;
            
        private:
            std::array<std::pair<NumberType,NumberType>,dim> _lineXD; 
    };
}
