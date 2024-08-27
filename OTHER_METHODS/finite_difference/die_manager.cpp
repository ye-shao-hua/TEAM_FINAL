
/**
 * 目前已经完成金属部分读取，接下来需要完成的是对于die的读取
 * polygon的vector已经初始化完毕，需要在txtread类中添加push_back
 * 方法用于将 读取die函数 中的polygon（还未写）添加到vector中
 * 读取完毕后还需写一个遍历函数循环整个vector，所以txtread类
 * 还需写一个begin与end方便遍历。
 *
 * 所需函数：
 *   diereader
 *   push_back
 *   begin
 *   end
 *   for_polygon
 *
 **/

#include <cstddef>
#include <die_manager/die_manager.h>
#include <fstream>
#include <iostream>
#include <optional>
#include <string_view>
#include <vector>

/**
 * 类 class AxxxBxxxCxxx
 * 函数 axxxBxxxCxxx
 * 私有成员变量_axxxx_bxxx_cxxx
 * 共有成员变量axxx_bxxx_cxxx
 * 计数变量 i j k n
 * 临时变量 buffer_xxx tmp_xxx
 * 枚举 常量 AXXX_BXXX_CXXX
 * 命名空间 随意
 */

// 命名空间
namespace bg = boost::geometry;
using point_2d = bg::model::d2::point_xy<double>;
using polygon_2d = bg::model::polygon<point_2d>;

void DieManager::open(std::string_view filename) {
  // 打开文件
  _ifs.open(filename.data());
  if (_ifs.is_open()) {
    _read_until_die();
    _read_die();
  }
}

bool DieManager::is_open() const { return _ifs.is_open(); }

void DieManager::_read_until_die() {
  // 读取前面不需要的金属部分
  std::string buffer;
  std::size_t number = 0;
  _ifs >> buffer;
  _ifs >> number;
  for (auto i = 0; i < 2 * number; i++) {
    _ifs >> buffer;
  }
  _ifs >> number;
  for (auto i = number; i > 0; i--) {
    _ifs >> buffer;
    std::size_t number = 0;
    _ifs >> number;
    for (auto j = 2 * number; j > 0; j--) {
      _ifs >> buffer;
    }
  }
}

void DieManager::_read_die() {
  size_t number{0};
  _ifs >> number;
  for (int i = 0; i < number; ++i) {
    std::string name{};
    _ifs >> name;
    size_t n_point{0};
    _ifs >> n_point;
    if (n_point == 0) {
      std::cerr << "i: " << i << "\t读取多边形的点数为0\n";
    }
    polygon_2d geometry;
    double x{}, y{};
    _ifs >> x >> y;
    bg::append(geometry, point_2d{x, y});
    for (int j = 1; j < n_point; ++j) {
      double x{}, y{};
      _ifs >> x >> y;
      bg::append(geometry, point_2d{x, y});
    }
    bg::append(geometry, point_2d{x, y});
    er_t er{0};
    _ifs >> er;
    _poly_vector.emplace_back(geometry, er);
  }
}

const std::optional<DieManager::er_t>
DieManager::search_die(const point_2d &point) const {
  static std::vector<std::pair<polygon_2d, er_t>>::const_iterator cache{};

  // 缓存命中
  if (cache != decltype(cache){} && bg::within(point, cache->first)) {
    return cache->second;
  }

  const auto &pos = std::find_if(
      _poly_vector.begin(), _poly_vector.end(),
      [&point](const auto &it) { return bg::within(point, it.first); });
  if (_poly_vector.end() == pos)
    return {};

  // 缓存
  cache = pos;

  return pos->second;
}
