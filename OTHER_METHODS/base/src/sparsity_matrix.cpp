#include "../include/fem/sparsity_matrix.h"

namespace fem{
    SparsityMatrix::SparsityMatrix(std::size_t rows,std::size_t cols):_rows(rows),_cols(cols){}
    std::size_t SparsityMatrix::rows() const{return _rows;}
    std::size_t SparsityMatrix::cols() const{return _cols;}
    double &SparsityMatrix::operator()(std::size_t row,std::size_t col){
        std::size_t index=row*_cols+col;
        return _map[std::size_t(index)];
    }
    const double& SparsityMatrix::operator()(std::size_t row,std::size_t col)const{
        static double zero = 0;
        std::size_t index=row*_cols+col;
        auto place=_map.find(index);
        if(place==_map.end())
            return zero;
        return place->second;
    }
    void SparsityMatrix::print(std::ostream &os)const{
        for(auto i=0;i<_rows;++i){
            for(auto j=0;j<_cols;++j){
                os<<operator()(i,j)<<",";
            }
            os<<"\n";
        }
    }
    SparsityMatrix SparsityMatrix::operator+(const SparsityMatrix &that) const{
        SparsityMatrix result=*this;
        for(const auto &i:that._map){
            result._map[std::size_t(i.first)]+=i.second;
        }
        return result;
    }


    SparsityMatrix SparsityMatrix::operator*(const SparsityMatrix &that) const{
        SparsityMatrix result{rows(),that.cols()};
        for(auto j=0;j<that.cols();++j){
            double sum{0};
            int row{0}; 
            int row_now{0};
            int col_now{0};
            for(const auto &i:_map){
                row_now = i.first / cols();
                col_now = i.first % cols();
                if(row_now != row){
                    result._map[std::size_t(row*that.cols()+j)]=sum;
                    sum=0;
                    row=row_now;
                }
                    sum+=i.second*that(col_now,j);
            }
            result._map[std::size_t(row*that.cols()+j)]=sum;

        }
        return result;
    }



}
