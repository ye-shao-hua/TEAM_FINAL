#pragma once
#include <array>

namespace fem{

template<int dim,typename NumberType> class Point{
    public:
        Point(const std::array<NumberType,dim> &data);
        NumberType &get(int index);
        const NumberType &get(int index) const;

    private:
        std::array<NumberType,dim> _data{};
};

}
