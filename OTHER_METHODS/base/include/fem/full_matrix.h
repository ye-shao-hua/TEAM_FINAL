#pragma once
#include <ostream>
#include <vector>

namespace fem {
class FullMatrix {
public:
  FullMatrix(int rows, int cols);
  int rows() const;
  int cols() const;
  double &operator()(int row, int col);
  const double &operator()(int row, int col) const;
  void print(std::ostream &os) const;
  FullMatrix operator+(const FullMatrix &that) const;
  FullMatrix operator+(double number) const;
  FullMatrix operator-() const;
  FullMatrix operator-(const FullMatrix &that) const;
  FullMatrix operator-(double number) const;
  FullMatrix operator*(const FullMatrix &that) const;
  FullMatrix operator*(double number) const;
  FullMatrix operator~() const;

protected:
  std::vector<std::vector<double>> _data;
};
} // namespace fem
