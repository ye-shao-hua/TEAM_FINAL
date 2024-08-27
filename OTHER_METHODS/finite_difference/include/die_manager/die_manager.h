#include <boost/geometry.hpp>
#include <fstream>
#include <optional>
#include <string_view>
#include <utility>
#include <vector>

class DieManager {
  using er_t = float;
  using point_2d = boost::geometry::model::d2::point_xy<double>;
  using polygon_2d = boost::geometry::model::polygon<point_2d>;

public:
  void open(std::string_view);
  bool is_open() const;
  const std::optional<er_t> search_die(const point_2d &) const;

private:
  void _read_until_die();
  void _read_die();

private:
  std::ifstream _ifs;
  std::vector<std::pair<polygon_2d, er_t>> _poly_vector;
};
