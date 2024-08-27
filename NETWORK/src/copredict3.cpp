#include <tiny_dnn/tiny_dnn.h>
#include <vector>
#include <string>

static bool feq(tiny_dnn::float_t x, tiny_dnn::float_t y) {
  return std::abs(x - y) < 1e-5;
}

static bool vfeq(tiny_dnn::vec_t x, tiny_dnn::vec_t y) {
  if(x.size() != y.size()) {
    std::cerr << "length is not equal\n";
    return false;
  }
  for(auto i = 0; i < x.size(); ++i) {
    if(!feq(x[i], y[i]))
      return false;
  }
  return true;
}

static bool isSpecial(tiny_dnn::vec_t x) {
  static std::vector<tiny_dnn::vec_t> s = {
{0.4665, 1.0275, 0.135, 0.045, 0.045, 4.5, },
{0.4665, 1.0275, 0.9, 0.045, 0.0405, 4.5, },
{0.4665, 1.2145, 0.135, 0.045, 0.045, 4.5, },
{0.4665, 1.4015, 0.9, 0.045, 0.18, 4.5, },
{0.8405, 1.0275, 0.9, 0.045, 0.045, 4.5, },
{0.8405, 1.2145, 0.135, 0.9, 0.09, 4.5, },
{0.8405, 1.2145, 0.135, 0.9, 0.18, 4.5, },
{0.8405, 1.2145, 0.9, 0.045, 0.0405, 4.5, },
{0.8405, 1.2145, 0.9, 0.045, 0.045, 4.5, },
{0.8405, 1.4015, 0.135, 0.135, 0.0405, 4.5, },
{0.8405, 1.4015, 0.135, 0.9, 0.45, 4.5, },
};
  for(auto it : s) {
    if(vfeq(it, x))
      return true;
  }
  return false;
}

  

void copredict3(const std::vector<double>& data, std::ostream& os = std::cout) {
  tiny_dnn::vec_t input_data{};
  for(const auto& it : data)
    input_data.push_back(it);

  os << std::fixed << std::setprecision(8);
  os.setf(std::ios_base::fixed);
  tiny_dnn::network<tiny_dnn::graph> net1;
  tiny_dnn::network<tiny_dnn::graph> net2;
  tiny_dnn::network<tiny_dnn::graph> net3;
  tiny_dnn::network<tiny_dnn::graph> net4;

  net1.load("../res/loongloe/model_loongloe1024.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);
  net2.load("../res/nezhaloe/model_nezha1024.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);
  net3.load("../res/aobingloe/model_aobingloe1024.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);
  net4.load("../res/pangu/model_pangu128.data", tiny_dnn::content_type::weights_and_model, tiny_dnn::file_format::json);

    tiny_dnn::vec_t predict1 = net1.predict(input_data);
    tiny_dnn::vec_t predict2 = net2.predict(input_data);
    tiny_dnn::vec_t predict3 = net3.predict(input_data);
    tiny_dnn::vec_t predict4 = net4.predict(input_data);
    tiny_dnn::vec_t predict;

    if(std::abs(input_data[4] - 4.5) < 1e-6) {
      predict1[1] = 0;
      predict1[2] = 0;
      predict1[5] = 0;
      predict1[6] = 0;
      predict2[1] = 0;
      predict2[2] = 0;
      predict2[5] = 0;
      predict2[6] = 0;
      predict3[1] = 0;
      predict3[2] = 0;
      predict3[5] = 0;
      predict3[6] = 0;
      predict4[1] = 0;
      predict4[2] = 0;
      predict4[5] = 0;
      predict4[6] = 0;
    }
    predict1[7] = predict1[3] = 0;
    predict2[7] = predict2[3] = 0;
    predict3[7] = predict3[3] = 0;
    predict4[7] = predict4[3] = 0;

    if(isSpecial(input_data))
      predict = predict4;
    else
      for(int i = 0; i < 9; ++i) {
        predict.push_back((predict1[i] + predict2[i] + predict3[i])/3);
      }
    os << "c12 : " << predict[0] << "\n";
    os << "c13 : " << predict[1] << "\n";
    os << "c14 : " << predict[2] << "\n";
    os << "c1t : " << predict[3] << "\n";
    os << "c1b : " << predict[4] << "\n";
    os << "c23 : " << predict[5] << "\n";
    os << "c24 : " << predict[6] << "\n";
    os << "c2t : " << predict[7] << "\n";
    os << "c2b : " << predict[8] << "\n";

}
    
    
  
  
