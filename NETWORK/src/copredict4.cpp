#include <tiny_dnn/tiny_dnn.h>
#include <vector>
#include <string>

void copredict4(const std::vector<double>& input, std::ostream& os = std::cout) {
    
  tiny_dnn::vec_t input_data;
  std::vector<double> mean{-0.0051,-0.0056,0.1414,-0.0714,-0.0437,-0.0030,-0.0090,-0.0030,-7.5194e-4};
  std::vector<double> sd{0.0025,0.0030,0.0458,0.0306,0.0338,0.0015,0.0146,0.0062,0.0012};
  for(const auto& it : input)
    input_data.push_back(it);
  os << std::fixed << std::setprecision(8);
  os.setf(std::ios_base::fixed);
  tiny_dnn::network<tiny_dnn::graph> net;

  net.load("../res/newdata/diag/hanaloe16384.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);


    tiny_dnn::vec_t predict = net.predict(input_data);
    for(auto i=0;i<9;++i){
        predict[i]=predict[i]*sd[i]+mean[i];
    }
    os << "botleft : " << predict[0] << "\n";
    os << "botright : " << predict[1] << "\n";
    os << "c1 : " << predict[2] << "\n";
    os << "c1Env : " << predict[3] << "\n";
    os << "c2 : " << predict[4] << "\n";
    os << "c2Env : " << predict[5] << "\n";
    os << "d1 : " << predict[6] << "\n";
    os << "d2 : " << predict[7] << "\n";
    os << "d2Env : " << predict[8] << "\n";
}
    
    
  
  
