#include "../include/fem/sparse_matrix.h"
#include <algorithm>
#include <cmath>

namespace fem::ex {

SparseMatrix::SparseMatrix(int rows, int cols)
    : _row_manager(rows), _cols(cols) {}

std::size_t SparseMatrix::rows() const { return _row_manager.size(); }

std::size_t SparseMatrix::cols() const { return _cols; }

void SparseMatrix::add(int row, int col, double value) {
  if (std::abs(value) < 1e-10)
    return;
  auto pos = std::find_if(_row_manager[row].begin(), _row_manager[row].end(),
                          [&](const std::pair<index_type, index_type> &p) {
                            return p.first == col;
                          });
  if (pos != _row_manager[row].end()) {
    _data[pos->second] += value;
    return;
  }

  _row_manager[row].emplace_back(col, _data.size());
  _data.push_back(value);
}

double SparseMatrix::operator()(int row, int col) const {
  auto pos = std::find_if(_row_manager[row].begin(), _row_manager[row].end(),
                          [&](const std::pair<index_type, index_type> &p) {
                            return p.first == col;
                          });
  if (pos != _row_manager[row].end())
    return _data[pos->second];
  return 0;
}

SparseMatrix SparseMatrix::operator*(const SparseMatrix &that) const {
  SparseMatrix result(rows(), that.cols());
  for (std::size_t j = 0; j < that.cols(); ++j) {
    for (std::size_t i = 0; i < rows(); ++i) {
      double value = 0;
      for (auto &it : _row_manager[i]) {
        auto pos =
            std::find_if(that._row_manager[it.first].begin(),
                         that._row_manager[it.first].end(),
                         [&](const std::pair<index_type, index_type> &p) {
                           return p.first == j;
                         });
        if (pos != that._row_manager[it.first].end()) {
          value += _data[it.second] * that._data[pos->second];
        }
      }
      result.add(i, j, value);
    }
  }
  return result;
}

} // namespace fem::ex
