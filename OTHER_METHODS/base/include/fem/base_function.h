#pragma once
#include <array>
#include <vector>

namespace fem {

template<int dim, typename NumberType> class Point;

template <typename NumberType> NumberType AttenuationFunction1D(double x);

template <typename NumberType> class BaseFunction1D {
public:
  NumberType operator()(double x);
  BaseFunction1D(int degree);
  NumberType grad(double x);

private:
  int _degree;
};

template <int dim, typename NumberType> class BaseFunctionGroupXD {
public:
  BaseFunctionGroupXD(int max_degree); // degree == 0, size == 1
  NumberType operator()(int id,const Point<dim,NumberType> &position);
  Point<dim, NumberType> grad(int id,const Point<dim,NumberType> &position);
private:
  std::array<std::vector<BaseFunction1D<NumberType>>, dim> _base_functions;

  //_base_functions[0]  x
  //_base_functions[0][4] x^4
  //_base_functions[1]  y
};

} // namespace fem
