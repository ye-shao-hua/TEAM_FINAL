#pragma once
#include <cmath>
#include "tiny_dnn/util/util.h"

namespace tiny_dnn {
class our_err {
 public:
  static float_t f(const vec_t &y, const vec_t &t) {
    assert(y.size() == t.size());
    float_t d{0.0};
    float_t d1{0.0};

    for (size_t i = 0; i < y.size(); ++i) {
      if(std::abs(t[i] - 0) < 1e-10) // real data is 0
        d += 0;
      else
        d += std::abs((y[i] - t[i]));
      d1 += std::abs(t[i]);
    }
    return d /= (d1 + 1e-10);
  }

  static vec_t df(const vec_t &y, const vec_t &t) {
    assert(y.size() == t.size());
    float_t factor{0.0};
    for (size_t i = 0; i < y.size(); ++i) factor += abs(t[i]);
    vec_t d(t.size());

    for (size_t i = 0; i < y.size(); ++i) {
      if(std::abs(t[i] - 0) < 1e-10)
        d[i] = 0;
      else
        d[i] = sig(y[i] - t[i]) / factor;
    }

    return d;
  }

 private:
  static float_t sig(float_t x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }
};

}  // namespace tiny_dnn
