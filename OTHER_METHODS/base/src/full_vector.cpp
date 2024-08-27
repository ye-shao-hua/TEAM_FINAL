#include "../include/fem/full_vector.h"
#include <cmath>

namespace fem {

FullVector::FullVector(int size) : FullMatrix(size, 1) {}

FullVector::FullVector(const FullMatrix &that) : FullMatrix(that) {
  if (that.cols() != 1)
    throw "that.cols() != 1";
}

double FullVector::operator*(const FullVector &that) const {
  double result = 0;
  for (int i = 0; i < _data.size(); ++i)
    result += operator()(i, 0) * that(i, 0);
  return result;
}

FullVector FullVector::operator*(double number) const {
  FullVector result(rows());
  for (int i = 0; i < rows(); ++i)
    result(i, 0) = operator()(i, 0) * number;
  return result;
}

double FullVector::norm() const {
  double result = 0;
  for (int i = 0; i < _data.size(); ++i)
    result += operator()(i, 0) * operator()(i, 0);
  result = std::sqrt(result);
  return result;
}

}
// namespace fem
