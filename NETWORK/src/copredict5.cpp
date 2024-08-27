#include <tiny_dnn/tiny_dnn.h>
#include <vector>
#include <string>

void copredict5(const std::vector<double>& input, std::ostream& os = std::cout) {
    
  tiny_dnn::vec_t input_data;
  std::vector<double> mean{-0.002434,-0.002290,0.119189,-0.050216,-0.008093,-0.056134};
  std::vector<double> sd{0.003975,0.003806,0.047172,0.034038,0.002953,0.034737};
  for(const auto& it : input)
    input_data.push_back(it);
  input_data.pop_back();
  os << std::fixed << std::setprecision(8);
  os.setf(std::ios_base::fixed);
  tiny_dnn::network<tiny_dnn::graph> net;

  net.load("../res/newdata/plate2l/hanaloe1024.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);


    tiny_dnn::vec_t predict = net.predict(input_data);
    for(auto i=0;i<6;++i){
        predict[i]=predict[i]*sd[i]+mean[i];
    }
    os << "botleft : " << predict[0] << "\n";
    os << "botright : " << predict[1] << "\n";
    os << "c1 : " << predict[2] << "\n";
    os << "c2 : " << predict[3] << "\n";
    os << "c2e : " << predict[4] << "\n";
    os << "c3 : " << predict[5] << "\n";
}
    
    
  
  
