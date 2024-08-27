#include <tiny_dnn/tiny_dnn.h>
#include <vector>
#include <string>

void copredict7(const std::vector<double>& input, std::ostream& os = std::cout) {
    
  tiny_dnn::vec_t input_data;
  std::vector<double> mean{-0.002974,-0.003153,0.146014,-0.071942,-0.043887,-0.002793,-0.0099929,-0.009939,-0.00133383};
  std::vector<double> sd{0.003635,0.004016,0.045421,0.031023,0.034431,0.001949,0.014634,0.012232,0.002701};
  for(const auto& it : input)
    input_data.push_back(it);
  os << std::fixed << std::setprecision(8);
  os.setf(std::ios_base::fixed);
  tiny_dnn::network<tiny_dnn::graph> net;

  net.load("../res/newdata/stack/leafloe16384.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);


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
    os << "d1L : " << predict[7] << "\n";
    os << "d2 : " << predict[8] << "\n";
}
    
    
  
  
