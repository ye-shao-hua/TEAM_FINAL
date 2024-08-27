#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

#define LINE 324

int main(int argc, char **argv){
    if(argc!=3){
        std::cout << "\n\nplease input table file (e.g. SUB-metal1_PLATE2L)\n\n"
            <<"and then input the output file name(e.g. a.txt) \n\n";
        return 0;
    }
    std::ifstream ifs{argv[1]};
    std::ofstream ofs{argv[2]};
    std::string buffer;
    double number;
    if(!ifs.is_open())
        std::cout<<"ifs open fail\n";
    if(!ofs.is_open())
        std::cout<<"ofs open fail\n";
    int number_of_line=0;
    int number_index=0;
get_line:
    getline(ifs,buffer);
    number_of_line=0;
    number_index=0;
    for(std::istringstream iss{buffer};!iss.eof();number_of_line++){
        iss>>buffer;
        if(buffer=="|")
            number_index=number_of_line;
    }
    if(number_index==0)
        goto get_line;
    std::cout<<"finished\n";
    std::cout<<number_of_line<<" "<<number_index<<"\n";
    ofs <<"index\t " <<"|\t " <<"c12\t " <<"c1t\t "<< "c1b\t "<<"c1tb\t "<<"c1e\t \n";
    for(auto i=0;i<LINE;++i){
        ofs<<i+1<<"\t "<< "|\t ";
        for(auto j=0;j<number_index+1;++j)
            ifs>>buffer;
        for(auto j=0;j<number_of_line-number_index-1;++j){
            ifs>>number;
            ofs<<number<<"\t ";
        }
        ofs<<"\n";
    }
    ifs.close();
    ofs.close();
    return 0;
}
