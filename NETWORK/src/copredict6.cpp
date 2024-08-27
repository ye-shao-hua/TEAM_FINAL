#include <tiny_dnn/tiny_dnn.h>
#include <vector>
#include <string>

void copredict6(const std::vector<double>& input, std::ostream& os = std::cout) {
    
  tiny_dnn::vec_t input_data;
  std::vector<double> mean{-0.001628,-0.001569,0.119803,-0.043488,-0.001686,-0.044645,-0.013572,-0.013209};
  std::vector<double> sd{0.003405,0.003310,0.044635,0.0350035,0.001553,0.035522,0.0047048,0.004936};
  for(const auto& it : input)
    input_data.push_back(it);
  input_data.pop_back();
  os << std::fixed << std::setprecision(8);
  os.setf(std::ios_base::fixed);
  tiny_dnn::network<tiny_dnn::graph> net;

  net.load("../res/newdata/plate3l/hanaloe2048.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);


    tiny_dnn::vec_t predict = net.predict(input_data);
    for(auto i=0;i<8;++i){
        predict[i]=predict[i]*sd[i]+mean[i];
    }
    os << "botleft : " << predict[0] << "\n";
    os << "botright : " << predict[1] << "\n";
    os << "c1 : " << predict[2] << "\n";
    os << "c2 : " << predict[3] << "\n";
    os << "c2e : " << predict[4] << "\n";
    os << "c3 : " << predict[5] << "\n";
    os << "topleft : " << predict[6] << "\n";
    os << "topright : " << predict[7] << "\n";
}
    
    
  
  
