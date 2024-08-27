#include <iostream>
#include <fstream>
#include <string>
int main(int argc, char **argv){
    if(argc!=4){
        std::cout<<"\n\nplease input out file name (e.g. out_1.txt)\nThen input number of readtimes(e.g. 5)\nLast input the index(e.g. 10)\n\n";
    return 0;
    }
    std::string buffer;
    double number;
    std::ifstream ifs{argv[1]};
    if(!ifs.is_open()){
        std::cout<<"file open fail\n";
        return 0;
    }
    std::cout<<std::stod(argv[3])<<"\t "<<"|\t ";
    for(auto i=0;i<std::stoi(argv[2]);++i){
        ifs>>buffer;
        ifs>>buffer;
        ifs>>number;
        std::cout<<number<<"\t ";
    }
    std::cout<<"\n";
    ifs.close();
    return 0;
}
