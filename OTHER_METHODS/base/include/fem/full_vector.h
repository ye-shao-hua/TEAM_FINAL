#pragma once
#include "./full_matrix.h"

namespace fem {

class FullVector : public FullMatrix {
public:
  explicit FullVector(int size);
  FullVector(const FullMatrix &that);
  double operator*(const FullVector &that) const;
  FullVector operator*(double number) const;
  double norm() const;
};

} // namespace fem
