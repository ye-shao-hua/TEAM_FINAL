#pragma once

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <tiny_dnn/layers/layer.h>
#include <tiny_dnn/util/util.h>

namespace tiny_dnn {

/**
 * element-wise log ```y = scale*log(x)+bias```
 **/
class log_layer : public layer {
 public:
  typedef layer Base;

  /**
   * @param in_shape [in] shape of input tensor
   * @param scale    [in] 
   * @param bias     [in] 
   */
  log_layer(const shape3d &in_shape,
            float_t scale = float_t{1.0},
            float_t bias = float_t{0.0})
    : layer({vector_type::data}, {vector_type::data}),
      in_shape_(in_shape),
      scale_(scale),
      bias_(bias) {}

  std::string layer_type() const override {return "log"; }

  std::vector<shape3d> in_shape() const override { return {in_shape_}; }

  std::vector<shape3d> out_shape() const override {return {in_shape_}; }

  void forward_propagation(const std::vector<tensor_t *> &in_data,
                           std::vector<tensor_t *> &out_data) override {
    const tensor_t &x = *in_data[0];
    tensor_t &y       = *out_data[0];

    for (size_t i = 0; i < x.size(); i++) {
      std::transform(x[i].begin(), x[i].end(), y[i].begin(),
                     [=](float_t x) { return scale_ * std::log(x) + bias_; });
    }
  }

  void back_propagation(const std::vector<tensor_t *> &in_data,
                        const std::vector<tensor_t *> &out_data,
                        std::vector<tensor_t *> &out_grad,
                        std::vector<tensor_t *> &in_grad) override {
    tensor_t &dx       = *in_grad[0];
    const tensor_t &dy = *out_grad[0];
    const tensor_t &x  = *in_data[0];
    (void)out_data;

    for (size_t i = 0; i < x.size(); i++) {
      for (size_t j = 0; j < x[i].size(); j++) {
        // f(x) = scale * log(x) + bias
        // ->
        //   dx = dy * df(x)
        //      = dy * scale / x
        dx[i][j] = dy[i][j] * scale_ / (x[i][j] + 1e-10);
      }
    }
  }

  float_t scale() const { return scale_; }

  float_t bias() const { return bias_; }

  friend struct serialization_buddy;

 private:
  shape3d in_shape_;
  float_t scale_;
  float_t bias_;
};

}  // namespace tiny_dnn
