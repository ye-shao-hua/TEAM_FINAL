#include "../include/fem/full_matrix.h"
#include <iostream>

namespace fem {

FullMatrix::FullMatrix(int rows, int cols)
    : _data(rows, std::vector<double>(cols, 0)) {}

int FullMatrix::rows() const { return _data.size(); }

int FullMatrix::cols() const {
  if (_data.size() == 0)
    throw "_data.size() == 0";
  return _data[0].size();
}

double &FullMatrix::operator()(int row, int col) { return _data[row][col]; }

const double &FullMatrix::operator()(int row, int col) const {
  return _data[row][col];
}

void FullMatrix::print(std::ostream &os) const {
  for (auto &it1 : _data) {
    for (auto &it2 : it1)
      os << it2 << ",";
    os << "\n";
  }
}

FullMatrix FullMatrix::operator+(const FullMatrix &that) const {
  FullMatrix result(_data.size(), _data[0].size());
  for (int i = 0; i < _data.size(); ++i)
    for (int j = 0; j < _data[0].size(); ++j)
      result(i, j) = operator()(i, j) + that(i, j);
  return result;
}

FullMatrix FullMatrix::operator+(double number) const {
  FullMatrix result(cols(), rows());
  for (int i = 0; i < cols(); ++i)
    for (int j = 0; j < rows(); ++j)
      result(i, j) = operator()(i, j) + number;
  return result;
}

FullMatrix FullMatrix::operator-() const {
  FullMatrix result(_data.size(), _data[0].size());
  for (int i = 0; i < _data.size(); ++i)
    for (int j = 0; j < _data[0].size(); ++j)
      result(i, j) = -operator()(i, j);
  return result;
}

FullMatrix FullMatrix::operator-(const FullMatrix &that) const {
  return operator+(-that);
}

FullMatrix FullMatrix::operator-(double number) const {
  return operator+(number);
}

FullMatrix FullMatrix::operator*(const FullMatrix &that) const {
  if (_data.size() == 0)
    throw "_data.size() == 0\n";
  if (that._data.size() == 0)
    throw "that._data.size() == 0\n";
  if (_data[0].size() != that._data.size())
    throw "this->_data.size() != that._data[0].size()\n";
  FullMatrix result(_data.size(), that._data[0].size());
  for (int i = 0; i < _data.size(); ++i) {           // row
    for (int j = 0; j < that._data[0].size(); ++j) { // col
      double sum = 0;
      for (int k = 0; k < _data[0].size(); ++k) {
        sum += _data[i][k] * that._data[k][j];
      }
      result(i, j) = sum;
    }
  }
  return result;
}

FullMatrix FullMatrix::operator*(double number) const {
  FullMatrix result(rows(), cols());
  for (int i = 0; i < cols(); ++i)
    for (int j = 0; j < rows(); ++j)
      result(i, j) = operator()(i, j) * number;
  return result;
}

FullMatrix FullMatrix::operator~() const {
  FullMatrix result(_data[0].size(), _data.size());
  for (int i = 0; i < _data.size(); ++i)
    for (int j = 0; j < _data[0].size(); ++j)
      result(j, i) = operator()(i, j);
  return result;
}

} // namespace fem
