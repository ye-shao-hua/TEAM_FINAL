#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
int main(int argc, char** argv){
    if(argc!=5){
        std::cout<<"\n\nerror! ! !\nPlease input realFileName, calculatedFileName,writeFileName and number of line\n\n";
        return 0;
    }
    std::ifstream real{argv[1]};
    std::ifstream calculated{argv[2]};
    std::ofstream write{argv[3]};
    int line{std::stoi(argv[4])};
    if(!real.is_open()){
        std::cout << "realFile open fail\n";
        return 0;
    }
    if(!calculated.is_open()){
        std::cout << "calculated file open fail\n";
        return 0;
    }
    if(!write.is_open()){
        std::cout << "writeFile open fail\n";
        return 0;
    }
    std::string buffer;
    std::string buffer_last;
    double number;
    getline(real,buffer);
    int index=0;
    int number_of_line=0;
    for(std::istringstream iss{buffer};!iss.eof();++number_of_line){
        iss>>buffer;
        if(buffer == "|")
            index=number_of_line;
        if(buffer==buffer_last){
            number_of_line--;
        }
        buffer_last=buffer;
    }
        std::cout <<index<<" "<<number_of_line<<"\n";
    int number_of_data=number_of_line-index-1;
    std::vector<double> error;
    double sum_real_cap=0;
    double sum_sub_cap=0;
    double real_buffer=0;
    double calculated_buffer;
    
    for(auto i=0;i<line;++i){
        real>>buffer;
        real>>buffer;
        calculated>>buffer;
        calculated>>buffer;
        sum_real_cap=0;
        sum_sub_cap=0;


        for(auto j=0;j<number_of_data;++j){
            real >> real_buffer;
            calculated >> calculated_buffer;
            sum_sub_cap+=std::fabs(calculated_buffer-real_buffer);
            sum_real_cap+=std::fabs(real_buffer);
        }
        error.push_back((sum_sub_cap/sum_real_cap)*100);
    }
	
    double average{0.};
    double max_value{0.};
    double min_value{999.};
    for(auto i:error){
        average+=i;
        if(i>max_value)
            max_value=i;
	if(i<min_value)
	    min_value=i;
    }
    average/=error.size();
    std::vector<double> sd(error);
    for(auto &i:sd){
	i=std::pow(i-average,2);	
    }
    double sum{0.};
    for(auto i:sd){
	sum+=i;	
    }
    sum /= error.size();
    sum=std::sqrt(sum);
    double two_sigma_down{0.};
    double two_sigma_up{0.};
    two_sigma_down = average-2*sum;
    two_sigma_up = average + 2* sum;
    
    write << "The average error:\t " << average << "%\t \n";
    write << "The max error:\t " << max_value << "%\t \n";
    write << "The min error:\t " << min_value << "%\t \n";
    write << "The standard error:\t " << sum << "%\t \n"; 
    write << "Two sigma region:{"<<two_sigma_down<<"%, "<<two_sigma_up << "%}\t \n";
    write << "\n";
    
    write << "every case's error:\t ";
    write << "\n";

    index = 1;
    for(auto i:error){
        write << index << " :\t ";
        write << "|\t ";
        write << i << "%\t \n";
        index++;
    }
    return 0;

}