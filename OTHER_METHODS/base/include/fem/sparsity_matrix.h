#pragma once
#include <map>
#include <ostream>

namespace fem {
    class SparsityMatrix{
        public:
            SparsityMatrix(std::size_t rows,std::size_t cols);
            std::size_t rows() const;
            std::size_t cols() const;
            double &operator()(std::size_t row, std::size_t col); // will add new element, if exist, recover. 
            const double& operator()(std::size_t row, std::size_t col) const; // read only, not add new element
            void print(std::ostream &os) const;
            SparsityMatrix operator+(const SparsityMatrix &that) const;
            SparsityMatrix operator*(const SparsityMatrix &that) const;

            private:
            std::size_t _rows,_cols;
            std::map<std::size_t,double> _map;
    };
}
