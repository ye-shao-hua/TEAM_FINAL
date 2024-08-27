#include <iostream>
#include <fstream>
#include <txtRead/txtRead.hpp>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Dense>
#include <ranges>
#include <chrono>

#define XLENGTH 1000  
//10000
#define YLENGTH 1000
#define K -16.5552  
//斜率

template<class FUNC>
void for_ysh(int begin,int end,FUNC func){
#pragma omp parallel for num_threads(8)
    for(auto i:std::ranges::iota_view{begin,end+1}){
        func(i);
    }
}

int main(){
    txtRead a;
    a.Open("b.txt");
    a.Read();
    //可直接操控单元

    //初始化数据
    double xmax=4.5,xmin=-4.5;
    double ymax=a.corner[3].y(),ymin=a.corner[0].y();
    Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(XLENGTH,YLENGTH); //初始化空间
    double dx=(xmax-xmin)/XLENGTH;
    double dy=(ymax-ymin)/YLENGTH;  //空间步长
    int xindex=0,yindex=0;//索引

    //映射关系
    //xindex=(0-xmin)/dx;
    //std::cout<<xindex<<" ";
    //std::cout<<xmin<<" "<<xmax<<" ";
    
    //函数声明
    void refresh(txtRead &a,Eigen::MatrixXf &mat);//刷新函数
    void calculate(txtRead& a,Eigen::MatrixXf &mat);   //计算函数
    
    //刷新计时测试
    /*
    for(auto i=0;i<10;i++){ 
        auto start=std::chrono::high_resolution_clock::now(); refresh(a,mat);
        auto end=std::chrono::high_resolution_clock::now();
        auto duration=std::chrono::duration_cast<std::chrono::nanoseconds>(end-start);
        std::cout<<duration.count()<<"\n";
    }*/

    //迭代计时测试
    /*
    refresh(a,mat);
    for(auto i=0;i<10;i++){ 
        auto start=std::chrono::high_resolution_clock::now();
        calculate(a,mat);
        auto end=std::chrono::high_resolution_clock::now();
        auto duration=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
        std::cout<<duration.count()<<"\n";
    }*/ 
    
    refresh(a,mat);

    auto start=std::chrono::high_resolution_clock::now();
    for(auto i=0;i<25000;i++){
        calculate(a,mat);
        refresh(a,mat);
    }
    auto end=std::chrono::high_resolution_clock::now();
    auto duration=std::chrono::duration_cast<std::chrono::seconds>(end-start);
    std::cout<<duration.count()<<"\n";



    //测试刷新后结果
    
    std::cout<<"正在保存"<<"\n";
    std::ofstream ofs{"out.csv"};
    for(auto i:std::ranges::iota_view{0,XLENGTH}){
        for(auto j:std::ranges::iota_view{0,YLENGTH}){
            ofs<<mat(i,j)<<",";
        }
        ofs<<"\b\n";
    }    
    

}

void refresh(txtRead &a,Eigen::MatrixXf& mat){//pattern3需要对分裂的金属特殊处理？or not
    double xmax=4.5,xmin=-4.5;
    double ymax=a.corner[3].y(),ymin=a.corner[0].y();
    double dx=(xmax-xmin)/XLENGTH;
    double dy=(ymax-ymin)/YLENGTH;  //空间步长
    //master
    double max=0;
    for(auto &i:a.master){
        if(i.y()>max){
            max=i.y();
        }
    }

    int max_i_index = static_cast<int>(std::ceil((max - a.master[0].y()) / dy));
    #pragma omp parallel for num_threads(8)
    for(int i_index = 0; i_index < max_i_index; ++i_index) {
        double i = a.master[0].y() + i_index * dy;
    //for(auto i=a.master[0].y();i<max;i=i+dy){
        int max_j_index = static_cast<int>(std::ceil((a.master[1].x()-(dy/K)*i_index - a.master[0].x() - (dy/K)*i_index) / dx));
        //#pragma omp parallel for num_threads(2)
        for(int j_index = 0; j_index < max_j_index; ++j_index) {
            double j = a.master[0].x() + (dy/K)*i_index + j_index * dx;
        //for(double j=a.master[0].x() + (dy/K)*n; j<a.master[1].x()-(dy/K)*n; j=j+dx){
            int xindex=(j-xmin)/dx;
            int yindex=(i-ymin)/dy;
            mat(xindex,yindex)=100;  //主导体电位恒为1
        }
    }
    //mental
    int n=0;
    #pragma omp parallel for num_threads(8)
    for(auto i:a.mental){
        if(i.get_name()!="botleft"&&i.get_name()!="botright"){
            max=0;
            for(auto j:i){
                if(j.y()>max)
                    max=j.y();
            }
            int max_j_index = static_cast<int>(std::ceil((max - i[0].y())/dy));
            for(int j_index = 0; j_index < max_j_index; ++j_index) {
                double j = i[0].y() + j_index * dy;
            // for(auto j=i[0].y();j<max;j=j+dy){
                for(auto k=i[0].x()+(dy/K)*n;k<i[1].x()-(dy/K)*n;k=k+dx){
                    int xindex=(k-xmin)/dx;
                    int yindex=(j-ymin)/dy;
                    mat(xindex,yindex)=10;  //非底部导体电位恒为0
                }
            n++;
            }
        }
    }
    #pragma omp parallel for num_threads(8)
    for(auto i=0;i<XLENGTH;i++){
    //for(auto i=xmin;i<xmax;i=i+dx){
        //int xindex=(i-xmin)/dx;
        //mat(xindex,1)=0;  //最下一行即可
        mat(i,0)=10; //底部导体电位恒为0
    }
    
}

void calculate(txtRead &a,Eigen::MatrixXf& mat){
    //先遍历四周，然后遍历内部  若矩阵空间为笛卡尔系，原点为矩阵左下，横轴为x，纵轴为y。
    //tip:下由于恒为1，不用给
    
    //计算值
    double xmax=4.5,xmin=-4.5;
    double ymax=a.corner[3].y(),ymin=a.corner[0].y();
    double dx=(xmax-xmin)/XLENGTH;
    double dy=(ymax-ymin)/YLENGTH;  //空间步长
    
    //上&左&右,边界条件全为零
    //上
    for_ysh(0,XLENGTH-1,[&mat](int i){
            mat(i,YLENGTH-1)=0;
    });
    //左&右
    for_ysh(1,YLENGTH-1,[&mat](int i){
        mat(1,i)=0;
        mat(XLENGTH-1,i)=0;
    });

    //中间
    for(auto i=1;i<XLENGTH-1;i++){
        for(auto j=1;j<YLENGTH-1;j++){
            mat(i,j)=(dy*dy*(mat(i+1,j)+mat(i-1,j))+dx*dx*(mat(i,j+1)+mat(i,j-1)))/(2*(dx*dx+dy*dy));
        }
    }
    /*
    for_ysh(1,XLENGTH-2,[&](int i){
        for_ysh(1,YLENGTH-2,[&](int j){
            mat(i,j)=(dy*dy*(mat(i+1,j)+mat(i-1,j))+dx*dx*(mat(i,j+1)+mat(i,j-1)))/(2*(dx*dx+dy*dy));
        });
    });*/
    
}





