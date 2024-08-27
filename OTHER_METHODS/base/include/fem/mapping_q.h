#pragma once
#include "full_matrix.h"
#include <vector>
#include "point.h"

namespace fem {

class FullMatrix;

template <int dim, typename NumberType> class MappingQ {
public:
  MappingQ(int degree);
  void set_jaccobi(int id, const FullMatrix &matrix);
  Point<dim,NumberType> mapping(Point<dim,NumberType> point);

private:
  std::vector<FullMatrix> _data;
};
} // namespace fem
