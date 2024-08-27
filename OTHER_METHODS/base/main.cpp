#include <iostream>
#include "./include/fem/base_function.h"
#include "include/fem/point.h"
#include "include/fem/sparsity_matrix.h"
#include "include/fem/sparse_matrix.h"
#include <chrono>

//#include "include/fem/full_matrix.h"
//#include "include/fem/mapping_q.h"
#include "include/fem/triangulation.h"
#include "include/fem/fe_q.h"
#include <array>
#include "include/fem/point.h"
#include "include/fem/cell.h"
#include "include/fem/hyper_cube.h"
using namespace fem;
int main(){
    /*
    FullMatrix fulla{1,1},fullb{1,1};
    fulla(0,0)=2;
    fullb(0,0)=4;
    MappingQ<1,double> mapa(2);
    mapa.set_jaccobi(0,fulla);
    mapa.set_jaccobi(1,fullb);

     
    Point<1,double> point1=mapa.mapping(Point<1,double>{std::array<double,1>{2.}});
    std::cout<<point1.get(0)<<" ";
    
    return 0;
    */
    FullMatrix full1(2,3),full2(2,2),full3(2,1);
    full1(0,0)=0;
    full1(0,1)=1;
    full1(0,2)=0;
    full1(1,0)=0;
    full1(1,1)=0;
    full1(1,2)=0;
    
    full2(0,0)=1;
    full2(0,1)=-1./2;
    full2(1,0)=0;
    full2(1,1)=1;

    full3(0,0)=0;
    full3(1,0)=0;


    full1.print(std::cout);
    std::cout<<"\n";
    full2.print(std::cout);
    std::cout<<"\n";
    full3.print(std::cout);
    std::cout<<"\n\n";

    MappingQ<2,double> map1(3);
    map1.set_jaccobi(0,full3);
    map1.set_jaccobi(1,full2);
    map1.set_jaccobi(2,full1);
    //std::cout<<full(1,1);

    Point<2,double> point=map1.mapping(Point<2,double>{std::array<double,2>{0.,0.}});
    std::cout<<"---------------------\n";
    std::cout<<point.get(0)<<" "<<point.get(1)<<"\n";
    return 0;

    FE_Q<2,double> fe_q{2};
    fe_q.init(1);
    Cell<2,double> ce{HyperCube<2, double>{{0.1,0.2},{0.1,0.2}}};
    fe_q.distribute(ce);

    std::cout<< ce(3,Point<2,double>{std::array<double,2>{0.1,0.2}}) <<"\n";
    return 0;


    auto start = std::chrono::steady_clock::now();
    SparsityMatrix sp1{100000,100000};
    SparsityMatrix sp12{100000,1};
    for(auto i=0;i<100000;++i)
        sp1(i, i) = 1.1;
    for(auto i=0;i<100000;++i)
        sp12(i, 0) = i;
    
    sp1=sp1*sp12;
    std::cout << sp1.rows() << " " << sp1.cols() << "\n";
    for(int i = 0; i < 10; ++i)
        std::cout << sp1(i, 0) << "\n";


    auto middle = std::chrono::steady_clock::now();
    ex::SparseMatrix sp2{100000,100000};
    for(auto i=0;i<100000;++i)
        sp2.add(i, i, 1.1);
    ex::SparseMatrix sp3{100000,1};
    for(auto i=0;i<100000;++i)
        sp3.add(i, 0, i);
    
    sp2=sp2*sp3;
    for(int i = 0; i < 10; ++i)
        std::cout << sp2(i, 0) << "\n";
    auto end = std::chrono::steady_clock::now();

    std::cout << (middle - start).count() << "\n";
    std::cout << (end - middle).count() << "\n";

    
    return 0;


    BaseFunction1D<double> a{4};
    BaseFunctionGroupXD<2, double> bf(2);
    Point<2, double> p = std::array<double ,2>{100., 5.};
    for(int i = 0; i < 8; ++i){
        for(int j = 0; j < 2; ++j) {
            std::cout << bf.grad(i, p).get(j) << " ";
        }
        std::cout << "\n";
    }
   // Point<2, double> grad = bf.grad(i, p);
    // 00 01 02
    // 10 11 12
    // 20 21 22
    // f = y * x^2
    //BaseFunctionGroupXD<1, typename NumberType>

    return 0;
}
