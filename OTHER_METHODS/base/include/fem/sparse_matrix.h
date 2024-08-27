#pragma once

#include <vector>
#include <map>

namespace fem::ex {

class SparseMatrix {
public:
  using index_type = int;

public:
  SparseMatrix(int rows, int cols);
  std::size_t rows() const;
  std::size_t cols() const;
  void add(int row, int col, double value);
  double operator()(int row, int col) const;
  SparseMatrix operator*(const SparseMatrix& that) const;

private:
  std::vector<double> _data;
  std::vector<std::vector<std::pair<index_type, index_type>>> _row_manager;
  std::size_t _cols;
};
} // namespace fem::ex
