#include "include/die_manager/die_manager.h"
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <fstream>
#include <iostream>
#include <ranges>
#include <txtRead/txtRead.hpp>
#include <solver/for_ysh.hpp>
#ifdef JC_USE_SSE
#include <xmmintrin.h> // 包含SSE指令集的头文件
#endif
#ifdef JC_USE_AVX
#undef JC_USE_SSE
#include <immintrin.h> // AVX
#endif

#define XLENGTH 1000
// 10000
#define YLENGTH 1000
// #define K -10.4803   //第一层
// #define K -16.5552     //其他层
//  斜率

// 用于选择横纵坐标范围
#define XMAX a.corner[1].x()
#define XMIN a.corner[0].x()
#define YMAX a.corner[3].y() / 2
#define YMIN a.corner[0].y()

/*
#define XMAX 4.5
#define XMIN -4.5
#define YMAX a.corner[3].y()
#define YMIN a.corner[0].y()
*/

void program_bar(float program) {
  std::cout << "\r";
  std::cout << "[";
  for (int i = 0; i < program * 100; ++i)
    std::cout << "#";
  for (int i = program * 100; i < 100; ++i)
    std::cout << " ";
  std::cout << "]" << std::flush;
}

int main() {
  DieManager d;
  d.open("a.txt");
  // std::cout << *d.search_die({0.5, 0.5}) << "\n";
  // return 0;

  txtRead a;
  a.Open("a.txt");
  a.Read();
  // 可直接操控单元

  // 初始化数据
  Eigen::MatrixXf mat = Eigen::MatrixXf::Zero(XLENGTH, YLENGTH); // 初始化空间

  // 映射关系
  // xindex=(0-xmin)/dx;
  // std::cout<<xindex<<" ";
  // std::cout<<xmin<<" "<<xmax<<" ";

  // 函数声明
  void refresh(txtRead & a, Eigen::MatrixXf & mat, double K = -10.); // 刷新函数
  void calculate(txtRead & a, Eigen::MatrixXf & mat, int iters = 1,
                 double K = -10.); // 计算函数
  void calculate_cap(txtRead & a, Eigen::MatrixXf & mat, DieManager & die,
                     double K = -10.);

  // 刷新计时测试
  /*
  for(auto i=0;i<10;i++){
      auto start=std::chrono::high_resolution_clock::now(); refresh(a,mat);
      auto end=std::chrono::high_resolution_clock::now();
      auto
  duration=std::chrono::duration_cast<std::chrono::nanoseconds>(end-start);
      std::cout<<duration.count()<<"\n";
  }*/

  // 迭代计时测试
  /*
  refresh(a,mat);
  for(auto i=0;i<10;i++){
      auto start=std::chrono::high_resolution_clock::now();
      calculate(a,mat);
      auto end=std::chrono::high_resolution_clock::now();
      auto
  duration=std::chrono::duration_cast<std::chrono::milliseconds>(end-start);
      std::cout<<duration.count()<<"\n";
  }*/
  double xmax = -10;
  for (auto i : a.master)
    if (i.x() > xmax)
      xmax = i.x();
  double ymax = 0;
  for (auto i : a.master)
    if (i.y() > ymax)
      ymax = i.y();
  double K = -(ymax - a.master[1].y()) / (xmax - a.master[1].x());
  std::cout << "斜率：" << K << "\n";

  refresh(a, mat, K);
  auto start = std::chrono::high_resolution_clock::now();
  for (auto i = 0; i < 1; i++) {
    calculate(a, mat, 100000, K);
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << duration.count() << "ms\n";

  // 测试刷新后结果

  std::cout << "正在保存"
            << "\n";
  std::ofstream ofs{"out.csv"};
  for (auto i : std::ranges::iota_view{0, XLENGTH}) {
    for (auto j : std::ranges::iota_view{0, YLENGTH}) {
      ofs << mat(i, j) << ",";
    }
    ofs << "\b\n";
  }

  // 计算电容
  calculate_cap(a, mat, d, K);
}
int n = 0;

void refresh(txtRead &a, Eigen::MatrixXf &mat,
             double K = -10.) { // pattern3需要对分裂的金属特殊处理？or not
  double xmax = XMAX, xmin = XMIN;
  double ymax = YMAX, ymin = YMIN;
  double dx = (xmax - xmin) / XLENGTH;
  double dy = (ymax - ymin) / YLENGTH; // 空间步长
  // master
  double max = 0;
  for (auto &i : a.master) {
    if (i.y() > max) {
      max = i.y();
    }
  } // 现在max存有master中最大y

  int max_i_index = static_cast<int>(
      std::ceil((max - a.master[0].y()) / dy)); // master中y最大索引
#pragma omp parallel for num_threads(8)
  for (int i_index = 0; i_index < max_i_index; ++i_index) {
    double i = a.master[0].y() + i_index * dy; // 还原，i等于y高度
    // for(auto i=a.master[0].y();i<max;i=i+dy){
    int max_j_index =
        static_cast<int>(std::ceil((a.master[1].x() - (dy / K) * i_index -
                                    a.master[0].x() - (dy / K) * i_index) /
                                   dx));
    // #pragma omp parallel for num_threads(2)
    for (int j_index = 0; j_index < max_j_index; ++j_index) {
      double j = a.master[0].x() + (dy / K) * i_index + j_index * dx;
      // for(double j=a.master[0].x() + (dy/K)*n; j<a.master[1].x()-(dy/K)*n;
      // j=j+dx){
      int xindex = (j - xmin) / dx;
      int yindex = (i - ymin) / dy;
      mat(xindex, yindex) = 1; // 主导体电位恒为1
    }
  }
  // mental
#pragma omp parallel for num_threads(8)
  for (auto i : a.mental) {
    if (i.get_name() != "botleft" && i.get_name() != "botright") {
      max = 0;
      for (auto j : i) {
        if (j.y() > max)
          max = j.y(); // 获得最大y
      }
      int n = 0;
      int max_j_index = static_cast<int>(std::ceil((max - i[0].y()) / dy));
      for (int j_index = 0; j_index < max_j_index; ++j_index) {
        double j = i[0].y() + j_index * dy;
        // for(auto j=i[0].y();j<max;j=j+dy){
        for (auto k = i[0].x() + (dy / K) * n; k < i[1].x() - (dy / K) * n;
             k = k + dx) {
          int xindex = (k - xmin) / dx;
          int yindex = (j - ymin) / dy;
          mat(xindex, yindex) = 0; // 非底部导体电位恒为0
        }
        n++;
      }
    }
  }
#pragma omp parallel for num_threads(8)
  for (auto i = 0; i < XLENGTH; i++) {
    // for(auto i=xmin;i<xmax;i=i+dx){
    // int xindex=(i-xmin)/dx;
    // mat(xindex,1)=0;  //最下一行即可
    mat(i, 0) = 0; // 底部导体电位恒为0
  }
  // 上&左&右,边界条件全为零
  // 上
  for_ysh(0, XLENGTH - 1, [&mat](int i) { mat(i, YLENGTH - 1) = 0; });
  // 左&右
  for_ysh(1, YLENGTH - 1, [&mat](int i) {
    mat(0, i) = 0;
    mat(XLENGTH - 1, i) = 0;
  });
}

void calculate_sse(txtRead &a, Eigen::MatrixXf &mat, int iters = 1,
                   double K = -10.) {
  std::clog << "调用calculate_sse()\n";
  double xmax = XMAX, xmin = XMIN;
  double ymax = YMAX, ymin = YMIN;
  double dx = (xmax - xmin) / XLENGTH;
  double dy = (ymax - ymin) / YLENGTH; // 空间步长

  Eigen::MatrixXf buffer = mat;
  refresh(a, buffer, K);

  for (int i = 0; i < iters; ++i) {
    program_bar(i / (float)iters);
    { // mat -> buffer
#pragma omp parallel for num_threads(8)
      for (int i = 1; i < XLENGTH - 1 - 8; i += 6) {
        int j;
        for (j = 1; j < YLENGTH - 1 - 8;
             j += 6) { // 注意j的步长为4，因为每次处理4个数据
          // 加载数据到SSE寄存器
          __m128 r0l = _mm_loadu_ps(&mat(i - 1, j - 1));
          __m128 r0r = _mm_loadu_ps(&mat(i + 3, j - 1));

          __m128 r1l = _mm_loadu_ps(&mat(i - 1, j));
          __m128 r1r = _mm_loadu_ps(&mat(i + 3, j));

          __m128 r2l = _mm_loadu_ps(&mat(i - 1, j + 1));
          __m128 r2r = _mm_loadu_ps(&mat(i + 3, j + 1));

          __m128 r3l = _mm_loadu_ps(&mat(i - 1, j + 2));
          __m128 r3r = _mm_loadu_ps(&mat(i + 3, j + 2));

          __m128 r4l = _mm_loadu_ps(&mat(i - 1, j + 3));
          __m128 r4r = _mm_loadu_ps(&mat(i + 3, j + 3));

          __m128 r5l = _mm_loadu_ps(&mat(i - 1, j + 4));
          __m128 r5r = _mm_loadu_ps(&mat(i + 3, j + 4));

          __m128 r6l = _mm_loadu_ps(&mat(i - 1, j + 5));
          __m128 r6r = _mm_loadu_ps(&mat(i + 3, j + 5));

          __m128 r7l = _mm_loadu_ps(&mat(i - 1, j + 6));
          __m128 r7r = _mm_loadu_ps(&mat(i + 3, j + 6));

// 对加载的数据进行转置
#define _MM_TRANS(r0l, r0r, r1l, r1r, r2l, r2r, r3l, r3r, r4l, r4r, r5l, r5r,  \
                  r6l, r6r, r7l, r7r, tr0l, tr0r, tr1l, tr1r, tr2l, tr2r,      \
                  tr3l, tr3r, tr4l, tr4r, tr5l, tr5r, tr6l, tr6r, tr7l, tr7r)  \
  {                                                                            \
    __m128 t1 = _mm_unpacklo_ps(r0l, r2l);                                     \
    __m128 t2 = _mm_unpacklo_ps(r0r, r2r);                                     \
    __m128 t3 = _mm_unpackhi_ps(r0l, r2l);                                     \
    __m128 t4 = _mm_unpackhi_ps(r0r, r2r);                                     \
                                                                               \
    __m128 t5 = _mm_unpacklo_ps(r1l, r3l);                                     \
    __m128 t6 = _mm_unpacklo_ps(r1r, r3r);                                     \
    __m128 t7 = _mm_unpackhi_ps(r1l, r3l);                                     \
    __m128 t8 = _mm_unpackhi_ps(r1r, r3r);                                     \
                                                                               \
    __m128 t9 = _mm_unpacklo_ps(r4l, r6l);                                     \
    __m128 t10 = _mm_unpacklo_ps(r4r, r6r);                                    \
    __m128 t11 = _mm_unpackhi_ps(r4l, r6l);                                    \
    __m128 t12 = _mm_unpackhi_ps(r4r, r6r);                                    \
                                                                               \
    __m128 t13 = _mm_unpacklo_ps(r5l, r7l);                                    \
    __m128 t14 = _mm_unpacklo_ps(r5r, r7r);                                    \
    __m128 t15 = _mm_unpackhi_ps(r5l, r7l);                                    \
    __m128 t16 = _mm_unpackhi_ps(r5r, r7r);                                    \
                                                                               \
    tr0l = _mm_unpacklo_ps(t1, t5);                                            \
    tr0r = _mm_unpacklo_ps(t9, t13);                                           \
    tr1l = _mm_unpackhi_ps(t1, t5);                                            \
    tr1r = _mm_unpackhi_ps(t9, t13);                                           \
    tr2l = _mm_unpacklo_ps(t3, t7);                                            \
    tr2r = _mm_unpacklo_ps(t11, t15);                                          \
    tr3l = _mm_unpackhi_ps(t3, t7);                                            \
    tr3r = _mm_unpackhi_ps(t11, t15);                                          \
                                                                               \
    tr4l = _mm_unpacklo_ps(t2, t6);                                            \
    tr4r = _mm_unpacklo_ps(t10, t14);                                          \
    tr5l = _mm_unpackhi_ps(t2, t6);                                            \
    tr5r = _mm_unpackhi_ps(t10, t14);                                          \
    tr6l = _mm_unpacklo_ps(t4, t8);                                            \
    tr6r = _mm_unpacklo_ps(t12, t16);                                          \
    tr7l = _mm_unpackhi_ps(t4, t8);                                            \
    tr7r = _mm_unpackhi_ps(t12, t16);                                          \
  }

#define _MM_SHOW(resultr1l)                                                    \
  for (int i = 0; i < 4; ++i) {                                                \
    std::cout << ((float *)&resultr1l)[i] << ", ";                             \
  }

          __m128 tr0l, tr0r, tr1l, tr1r, tr2l, tr2r, tr3l, tr3r, tr4l, tr4r,
              tr5l, tr5r, tr6l, tr6r, tr7l, tr7r;
          _MM_TRANS(r0l, r0r, r1l, r1r, r2l, r2r, r3l, r3r, r4l, r4r, r5l, r5r,
                    r6l, r6r, r7l, r7r, tr0l, tr0r, tr1l, tr1r, tr2l, tr2r,
                    tr3l, tr3r, tr4l, tr4r, tr5l, tr5r, tr6l, tr6r, tr7l, tr7r);

          // 算法
          __m128 sumr1l = _mm_add_ps(r0l, r2l);
          __m128 sumr1r = _mm_add_ps(r0r, r2r);
          __m128 sumr2l = _mm_add_ps(r1l, r3l);
          __m128 sumr2r = _mm_add_ps(r1r, r3r);
          __m128 sumr3l = _mm_add_ps(r2l, r4l);
          __m128 sumr3r = _mm_add_ps(r2r, r4r);
          __m128 sumr4l = _mm_add_ps(r3l, r5l);
          __m128 sumr4r = _mm_add_ps(r3r, r5r);
          __m128 sumr5l = _mm_add_ps(r4l, r6l);
          __m128 sumr5r = _mm_add_ps(r4r, r6r);
          __m128 sumr6l = _mm_add_ps(r5l, r7l);
          __m128 sumr6r = _mm_add_ps(r5r, r7r);

          __m128 sumtr1l = _mm_add_ps(tr0l, tr2l);
          __m128 sumtr1r = _mm_add_ps(tr0r, tr2r);
          __m128 sumtr2l = _mm_add_ps(tr1l, tr3l);
          __m128 sumtr2r = _mm_add_ps(tr1r, tr3r);
          __m128 sumtr3l = _mm_add_ps(tr2l, tr4l);
          __m128 sumtr3r = _mm_add_ps(tr2r, tr4r);
          __m128 sumtr4l = _mm_add_ps(tr3l, tr5l);
          __m128 sumtr4r = _mm_add_ps(tr3r, tr5r);
          __m128 sumtr5l = _mm_add_ps(tr4l, tr6l);
          __m128 sumtr5r = _mm_add_ps(tr4r, tr6r);
          __m128 sumtr6l = _mm_add_ps(tr5l, tr7l);
          __m128 sumtr6r = _mm_add_ps(tr5r, tr7r);

          const __m128 dx2 = _mm_set1_ps(dx * dx);
          const __m128 dy2 = _mm_set1_ps(dy * dy);

          __m128 dx2_sumr1l = _mm_mul_ps(dx2, sumr1l);
          __m128 dx2_sumr1r = _mm_mul_ps(dx2, sumr1r);
          __m128 dx2_sumr2l = _mm_mul_ps(dx2, sumr2l);
          __m128 dx2_sumr2r = _mm_mul_ps(dx2, sumr2r);
          __m128 dx2_sumr3l = _mm_mul_ps(dx2, sumr3l);
          __m128 dx2_sumr3r = _mm_mul_ps(dx2, sumr3r);
          __m128 dx2_sumr4l = _mm_mul_ps(dx2, sumr4l);
          __m128 dx2_sumr4r = _mm_mul_ps(dx2, sumr4r);
          __m128 dx2_sumr5l = _mm_mul_ps(dx2, sumr5l);
          __m128 dx2_sumr5r = _mm_mul_ps(dx2, sumr5r);
          __m128 dx2_sumr6l = _mm_mul_ps(dx2, sumr6l);
          __m128 dx2_sumr6r = _mm_mul_ps(dx2, sumr6r);

          __m128 dy2_sumtr1l = _mm_mul_ps(dy2, sumtr1l);
          __m128 dy2_sumtr1r = _mm_mul_ps(dy2, sumtr1r);
          __m128 dy2_sumtr2l = _mm_mul_ps(dy2, sumtr2l);
          __m128 dy2_sumtr2r = _mm_mul_ps(dy2, sumtr2r);
          __m128 dy2_sumtr3l = _mm_mul_ps(dy2, sumtr3l);
          __m128 dy2_sumtr3r = _mm_mul_ps(dy2, sumtr3r);
          __m128 dy2_sumtr4l = _mm_mul_ps(dy2, sumtr4l);
          __m128 dy2_sumtr4r = _mm_mul_ps(dy2, sumtr4r);
          __m128 dy2_sumtr5l = _mm_mul_ps(dy2, sumtr5l);
          __m128 dy2_sumtr5r = _mm_mul_ps(dy2, sumtr5r);
          __m128 dy2_sumtr6l = _mm_mul_ps(dy2, sumtr6l);
          __m128 dy2_sumtr6r = _mm_mul_ps(dy2, sumtr6r);
          __m128 dy2_reverse7l;
          __m128 dy2_reverse7r;
          __m128 dy2_reverse8l;
          __m128 dy2_reverse8r;

          __m128 dy2_sumttr1l;
          __m128 dy2_sumttr1r;
          __m128 dy2_sumttr2l;
          __m128 dy2_sumttr2r;
          __m128 dy2_sumttr3l;
          __m128 dy2_sumttr3r;
          __m128 dy2_sumttr4l;
          __m128 dy2_sumttr4r;
          __m128 dy2_sumttr5l;
          __m128 dy2_sumttr5r;
          __m128 dy2_sumttr6l;
          __m128 dy2_sumttr6r;
          __m128 dy2_reverset7l;
          __m128 dy2_reverset7r;
          __m128 dy2_reverset8l;
          __m128 dy2_reverset8r;

          _MM_TRANS(dy2_reverse8l, dy2_reverse8r, dy2_sumtr1l, dy2_sumtr1r,
                    dy2_sumtr2l, dy2_sumtr2r, dy2_sumtr3l, dy2_sumtr3r,
                    dy2_sumtr4l, dy2_sumtr4r, dy2_sumtr5l, dy2_sumtr5r,
                    dy2_sumtr6l, dy2_sumtr6r, dy2_reverse7l, dy2_reverse7r,

                    dy2_reverset8l, dy2_reverset8r, dy2_sumttr1l, dy2_sumttr1r,
                    dy2_sumttr2l, dy2_sumttr2r, dy2_sumttr3l, dy2_sumttr3r,
                    dy2_sumttr4l, dy2_sumttr4r, dy2_sumttr5l, dy2_sumttr5r,
                    dy2_sumttr6l, dy2_sumttr6r, dy2_reverset7l, dy2_reverset7r);

          const __m128 $2dx2_dy2 =
              _mm_mul_ps(_mm_add_ps(dx2, dy2), _mm_set1_ps(2));
          __m128 resultr1l =
              _mm_div_ps(_mm_add_ps(dx2_sumr1l, dy2_sumttr1l), $2dx2_dy2);
          __m128 resultr1r =
              _mm_div_ps(_mm_add_ps(dx2_sumr1r, dy2_sumttr1r), $2dx2_dy2);
          __m128 resultr2l =
              _mm_div_ps(_mm_add_ps(dx2_sumr2l, dy2_sumttr2l), $2dx2_dy2);
          __m128 resultr2r =
              _mm_div_ps(_mm_add_ps(dx2_sumr2r, dy2_sumttr2r), $2dx2_dy2);
          __m128 resultr3l =
              _mm_div_ps(_mm_add_ps(dx2_sumr3l, dy2_sumttr3l), $2dx2_dy2);
          __m128 resultr3r =
              _mm_div_ps(_mm_add_ps(dx2_sumr3r, dy2_sumttr3r), $2dx2_dy2);
          __m128 resultr4l =
              _mm_div_ps(_mm_add_ps(dx2_sumr4l, dy2_sumttr4l), $2dx2_dy2);
          __m128 resultr4r =
              _mm_div_ps(_mm_add_ps(dx2_sumr4r, dy2_sumttr4r), $2dx2_dy2);
          __m128 resultr5l =
              _mm_div_ps(_mm_add_ps(dx2_sumr5l, dy2_sumttr5l), $2dx2_dy2);
          __m128 resultr5r =
              _mm_div_ps(_mm_add_ps(dx2_sumr5r, dy2_sumttr5r), $2dx2_dy2);
          __m128 resultr6l =
              _mm_div_ps(_mm_add_ps(dx2_sumr6l, dy2_sumttr6l), $2dx2_dy2);
          __m128 resultr6r =
              _mm_div_ps(_mm_add_ps(dx2_sumr6r, dy2_sumttr6r), $2dx2_dy2);

          __m128 result_reverse1l, result_reverse1r, result_reverse2l,
              result_reverse2r;
          __m128 t_resultr1l, t_resultr1r, t_resultr2l, t_resultr2r,
              t_resultr3l, t_resultr3r, t_resultr4l, t_resultr4r, t_resultr5l,
              t_resultr5r, t_resultr6l, t_resultr6r, t_resultr7l, t_resultr7r,
              t_resultr8l, t_resultr8r;
          _MM_TRANS(resultr1l, resultr1r, resultr2l, resultr2r, resultr3l,
                    resultr3r, resultr4l, resultr4r, resultr5l, resultr5r,
                    resultr6l, resultr6r, result_reverse1l, result_reverse1r,
                    result_reverse2l, result_reverse2r,

                    t_resultr1l, t_resultr1r, t_resultr2l, t_resultr2r,
                    t_resultr3l, t_resultr3r, t_resultr4l, t_resultr4r,
                    t_resultr5l, t_resultr5r, t_resultr6l, t_resultr6r,
                    t_resultr7l, t_resultr7r, t_resultr8l, t_resultr8r);
          _MM_TRANS(t_resultr2l, t_resultr2r, t_resultr3l, t_resultr3r,
                    t_resultr4l, t_resultr4r, t_resultr5l, t_resultr5r,
                    t_resultr6l, t_resultr6r, t_resultr7l, t_resultr7r,
                    t_resultr8l, t_resultr8r, result_reverse1l,
                    result_reverse1r,

                    resultr1l, resultr1r, resultr2l, resultr2r, resultr3l,
                    resultr3r, resultr4l, resultr4r, resultr5l, resultr5r,
                    resultr6l, resultr6r, result_reverse1l, result_reverse1r,
                    result_reverse2l, result_reverse2r);

          _mm_storeu_ps(&buffer(i, j), resultr1l);
          _mm_storel_epi64((__m128i_u *)&buffer(i + 4, j),
                           *(__m128i *)&resultr1r);
          //_mm_storeu_ps(&buffer(i+4, j), resultr1r);
          _mm_storeu_ps(&buffer(i, j + 1), resultr2l);
          _mm_storel_epi64((__m128i_u *)&buffer(i + 4, j + 1),
                           *(__m128i *)&resultr2r);
          //_mm_storeu_ps(&buffer(i+4, j+1), resultr2r);
          _mm_storeu_ps(&buffer(i, j + 2), resultr3l);
          _mm_storel_epi64((__m128i_u *)&buffer(i + 4, j + 2),
                           *(__m128i *)&resultr3r);
          //_mm_storeu_ps(&buffer(i+4, j+2), resultr3r);
          _mm_storeu_ps(&buffer(i, j + 3), resultr4l);
          _mm_storel_epi64((__m128i_u *)&buffer(i + 4, j + 3),
                           *(__m128i *)&resultr4r);
          //_mm_storeu_ps(&buffer(i+4, j+3), resultr4r);
          _mm_storeu_ps(&buffer(i, j + 4), resultr5l);
          _mm_storel_epi64((__m128i_u *)&buffer(i + 4, j + 4),
                           *(__m128i *)&resultr5r);
          //_mm_storeu_ps(&buffer(i+4, j+4), resultr5r);
          _mm_storeu_ps(&buffer(i, j + 5), resultr6l);
          _mm_storel_epi64((__m128i_u *)&buffer(i + 4, j + 5),
                           *(__m128i *)&resultr6r);
          //_mm_storeu_ps(&buffer(i+4, j+5), resultr6r);
        }
        for (; j < YLENGTH - 1; ++j) {
          for (int ii = i; ii < i + 6; ++ii) {
            buffer(ii, j) = (dy * dy * (mat(ii + 1, j) + mat(ii - 1, j)) +
                             dx * dx * (mat(ii, j + 1) + mat(ii, j - 1))) /
                            (2 * (dx * dx + dy * dy));
          }
        }
        for (; i + 6 >= XLENGTH - 1 - 8 && i < XLENGTH - 1; i += 1) {
          for (int j = 1; j < YLENGTH - 1; ++j) {
            buffer(i, j) = (dy * dy * (mat(i + 1, j) + mat(i - 1, j)) +
                            dx * dx * (mat(i, j + 1) + mat(i, j - 1))) /
                           (2 * (dx * dx + dy * dy));
          }
        }
      }
    } // end mat -> buffer
    refresh(a, buffer, K);

    { // buffer -> mat
#pragma omp parallel for num_threads(8)
      for (int i = 1; i < XLENGTH - 1 - 8; i += 6) {
        int j;
        for (j = 1; j < YLENGTH - 1 - 8;
             j += 6) { // 注意j的步长为4，因为每次处理4个数据
          // 加载数据到SSE寄存器
          __m128 r0l = _mm_loadu_ps(&buffer(i - 1, j - 1));
          __m128 r0r = _mm_loadu_ps(&buffer(i + 3, j - 1));

          __m128 r1l = _mm_loadu_ps(&buffer(i - 1, j));
          __m128 r1r = _mm_loadu_ps(&buffer(i + 3, j));

          __m128 r2l = _mm_loadu_ps(&buffer(i - 1, j + 1));
          __m128 r2r = _mm_loadu_ps(&buffer(i + 3, j + 1));

          __m128 r3l = _mm_loadu_ps(&buffer(i - 1, j + 2));
          __m128 r3r = _mm_loadu_ps(&buffer(i + 3, j + 2));

          __m128 r4l = _mm_loadu_ps(&buffer(i - 1, j + 3));
          __m128 r4r = _mm_loadu_ps(&buffer(i + 3, j + 3));

          __m128 r5l = _mm_loadu_ps(&buffer(i - 1, j + 4));
          __m128 r5r = _mm_loadu_ps(&buffer(i + 3, j + 4));

          __m128 r6l = _mm_loadu_ps(&buffer(i - 1, j + 5));
          __m128 r6r = _mm_loadu_ps(&buffer(i + 3, j + 5));

          __m128 r7l = _mm_loadu_ps(&buffer(i - 1, j + 6));
          __m128 r7r = _mm_loadu_ps(&buffer(i + 3, j + 6));

          // 对加载的数据进行转置

          __m128 tr0l, tr0r, tr1l, tr1r, tr2l, tr2r, tr3l, tr3r, tr4l, tr4r,
              tr5l, tr5r, tr6l, tr6r, tr7l, tr7r;
          _MM_TRANS(r0l, r0r, r1l, r1r, r2l, r2r, r3l, r3r, r4l, r4r, r5l, r5r,
                    r6l, r6r, r7l, r7r, tr0l, tr0r, tr1l, tr1r, tr2l, tr2r,
                    tr3l, tr3r, tr4l, tr4r, tr5l, tr5r, tr6l, tr6r, tr7l, tr7r);

          // 算法
          __m128 sumr1l = _mm_add_ps(r0l, r2l);
          __m128 sumr1r = _mm_add_ps(r0r, r2r);
          __m128 sumr2l = _mm_add_ps(r1l, r3l);
          __m128 sumr2r = _mm_add_ps(r1r, r3r);
          __m128 sumr3l = _mm_add_ps(r2l, r4l);
          __m128 sumr3r = _mm_add_ps(r2r, r4r);
          __m128 sumr4l = _mm_add_ps(r3l, r5l);
          __m128 sumr4r = _mm_add_ps(r3r, r5r);
          __m128 sumr5l = _mm_add_ps(r4l, r6l);
          __m128 sumr5r = _mm_add_ps(r4r, r6r);
          __m128 sumr6l = _mm_add_ps(r5l, r7l);
          __m128 sumr6r = _mm_add_ps(r5r, r7r);

          __m128 sumtr1l = _mm_add_ps(tr0l, tr2l);
          __m128 sumtr1r = _mm_add_ps(tr0r, tr2r);
          __m128 sumtr2l = _mm_add_ps(tr1l, tr3l);
          __m128 sumtr2r = _mm_add_ps(tr1r, tr3r);
          __m128 sumtr3l = _mm_add_ps(tr2l, tr4l);
          __m128 sumtr3r = _mm_add_ps(tr2r, tr4r);
          __m128 sumtr4l = _mm_add_ps(tr3l, tr5l);
          __m128 sumtr4r = _mm_add_ps(tr3r, tr5r);
          __m128 sumtr5l = _mm_add_ps(tr4l, tr6l);
          __m128 sumtr5r = _mm_add_ps(tr4r, tr6r);
          __m128 sumtr6l = _mm_add_ps(tr5l, tr7l);
          __m128 sumtr6r = _mm_add_ps(tr5r, tr7r);

          const __m128 dx2 = _mm_set1_ps(dx * dx);
          const __m128 dy2 = _mm_set1_ps(dy * dy);

          __m128 dx2_sumr1l = _mm_mul_ps(dx2, sumr1l);
          __m128 dx2_sumr1r = _mm_mul_ps(dx2, sumr1r);
          __m128 dx2_sumr2l = _mm_mul_ps(dx2, sumr2l);
          __m128 dx2_sumr2r = _mm_mul_ps(dx2, sumr2r);
          __m128 dx2_sumr3l = _mm_mul_ps(dx2, sumr3l);
          __m128 dx2_sumr3r = _mm_mul_ps(dx2, sumr3r);
          __m128 dx2_sumr4l = _mm_mul_ps(dx2, sumr4l);
          __m128 dx2_sumr4r = _mm_mul_ps(dx2, sumr4r);
          __m128 dx2_sumr5l = _mm_mul_ps(dx2, sumr5l);
          __m128 dx2_sumr5r = _mm_mul_ps(dx2, sumr5r);
          __m128 dx2_sumr6l = _mm_mul_ps(dx2, sumr6l);
          __m128 dx2_sumr6r = _mm_mul_ps(dx2, sumr6r);

          __m128 dy2_sumtr1l = _mm_mul_ps(dy2, sumtr1l);
          __m128 dy2_sumtr1r = _mm_mul_ps(dy2, sumtr1r);
          __m128 dy2_sumtr2l = _mm_mul_ps(dy2, sumtr2l);
          __m128 dy2_sumtr2r = _mm_mul_ps(dy2, sumtr2r);
          __m128 dy2_sumtr3l = _mm_mul_ps(dy2, sumtr3l);
          __m128 dy2_sumtr3r = _mm_mul_ps(dy2, sumtr3r);
          __m128 dy2_sumtr4l = _mm_mul_ps(dy2, sumtr4l);
          __m128 dy2_sumtr4r = _mm_mul_ps(dy2, sumtr4r);
          __m128 dy2_sumtr5l = _mm_mul_ps(dy2, sumtr5l);
          __m128 dy2_sumtr5r = _mm_mul_ps(dy2, sumtr5r);
          __m128 dy2_sumtr6l = _mm_mul_ps(dy2, sumtr6l);
          __m128 dy2_sumtr6r = _mm_mul_ps(dy2, sumtr6r);
          __m128 dy2_reverse7l;
          __m128 dy2_reverse7r;
          __m128 dy2_reverse8l;
          __m128 dy2_reverse8r;

          __m128 dy2_sumttr1l;
          __m128 dy2_sumttr1r;
          __m128 dy2_sumttr2l;
          __m128 dy2_sumttr2r;
          __m128 dy2_sumttr3l;
          __m128 dy2_sumttr3r;
          __m128 dy2_sumttr4l;
          __m128 dy2_sumttr4r;
          __m128 dy2_sumttr5l;
          __m128 dy2_sumttr5r;
          __m128 dy2_sumttr6l;
          __m128 dy2_sumttr6r;
          __m128 dy2_reverset7l;
          __m128 dy2_reverset7r;
          __m128 dy2_reverset8l;
          __m128 dy2_reverset8r;

          _MM_TRANS(dy2_reverse8l, dy2_reverse8r, dy2_sumtr1l, dy2_sumtr1r,
                    dy2_sumtr2l, dy2_sumtr2r, dy2_sumtr3l, dy2_sumtr3r,
                    dy2_sumtr4l, dy2_sumtr4r, dy2_sumtr5l, dy2_sumtr5r,
                    dy2_sumtr6l, dy2_sumtr6r, dy2_reverse7l, dy2_reverse7r,

                    dy2_reverset8l, dy2_reverset8r, dy2_sumttr1l, dy2_sumttr1r,
                    dy2_sumttr2l, dy2_sumttr2r, dy2_sumttr3l, dy2_sumttr3r,
                    dy2_sumttr4l, dy2_sumttr4r, dy2_sumttr5l, dy2_sumttr5r,
                    dy2_sumttr6l, dy2_sumttr6r, dy2_reverset7l, dy2_reverset7r);

          const __m128 $2dx2_dy2 =
              _mm_mul_ps(_mm_add_ps(dx2, dy2), _mm_set1_ps(2));
          __m128 resultr1l =
              _mm_div_ps(_mm_add_ps(dx2_sumr1l, dy2_sumttr1l), $2dx2_dy2);
          __m128 resultr1r =
              _mm_div_ps(_mm_add_ps(dx2_sumr1r, dy2_sumttr1r), $2dx2_dy2);
          __m128 resultr2l =
              _mm_div_ps(_mm_add_ps(dx2_sumr2l, dy2_sumttr2l), $2dx2_dy2);
          __m128 resultr2r =
              _mm_div_ps(_mm_add_ps(dx2_sumr2r, dy2_sumttr2r), $2dx2_dy2);
          __m128 resultr3l =
              _mm_div_ps(_mm_add_ps(dx2_sumr3l, dy2_sumttr3l), $2dx2_dy2);
          __m128 resultr3r =
              _mm_div_ps(_mm_add_ps(dx2_sumr3r, dy2_sumttr3r), $2dx2_dy2);
          __m128 resultr4l =
              _mm_div_ps(_mm_add_ps(dx2_sumr4l, dy2_sumttr4l), $2dx2_dy2);
          __m128 resultr4r =
              _mm_div_ps(_mm_add_ps(dx2_sumr4r, dy2_sumttr4r), $2dx2_dy2);
          __m128 resultr5l =
              _mm_div_ps(_mm_add_ps(dx2_sumr5l, dy2_sumttr5l), $2dx2_dy2);
          __m128 resultr5r =
              _mm_div_ps(_mm_add_ps(dx2_sumr5r, dy2_sumttr5r), $2dx2_dy2);
          __m128 resultr6l =
              _mm_div_ps(_mm_add_ps(dx2_sumr6l, dy2_sumttr6l), $2dx2_dy2);
          __m128 resultr6r =
              _mm_div_ps(_mm_add_ps(dx2_sumr6r, dy2_sumttr6r), $2dx2_dy2);

          __m128 result_reverse1l, result_reverse1r, result_reverse2l,
              result_reverse2r;
          __m128 t_resultr1l, t_resultr1r, t_resultr2l, t_resultr2r,
              t_resultr3l, t_resultr3r, t_resultr4l, t_resultr4r, t_resultr5l,
              t_resultr5r, t_resultr6l, t_resultr6r, t_resultr7l, t_resultr7r,
              t_resultr8l, t_resultr8r;
          _MM_TRANS(resultr1l, resultr1r, resultr2l, resultr2r, resultr3l,
                    resultr3r, resultr4l, resultr4r, resultr5l, resultr5r,
                    resultr6l, resultr6r, result_reverse1l, result_reverse1r,
                    result_reverse2l, result_reverse2r,

                    t_resultr1l, t_resultr1r, t_resultr2l, t_resultr2r,
                    t_resultr3l, t_resultr3r, t_resultr4l, t_resultr4r,
                    t_resultr5l, t_resultr5r, t_resultr6l, t_resultr6r,
                    t_resultr7l, t_resultr7r, t_resultr8l, t_resultr8r);

          _MM_TRANS(t_resultr2l, t_resultr2r, t_resultr3l, t_resultr3r,
                    t_resultr4l, t_resultr4r, t_resultr5l, t_resultr5r,
                    t_resultr6l, t_resultr6r, t_resultr7l, t_resultr7r,
                    t_resultr8l, t_resultr8r, result_reverse1l,
                    result_reverse1r,

                    resultr1l, resultr1r, resultr2l, resultr2r, resultr3l,
                    resultr3r, resultr4l, resultr4r, resultr5l, resultr5r,
                    resultr6l, resultr6r, result_reverse1l, result_reverse1r,
                    result_reverse2l, result_reverse2r);

          _mm_storeu_ps(&mat(i, j), resultr1l);
          _mm_storel_epi64((__m128i_u *)&mat(i + 4, j), *(__m128i *)&resultr1r);

          _mm_storeu_ps(&mat(i, j + 1), resultr2l);
          _mm_storel_epi64((__m128i_u *)&mat(i + 4, j + 1),
                           *(__m128i *)&resultr2r);
          _mm_storeu_ps(&mat(i, j + 2), resultr3l);
          _mm_storel_epi64((__m128i_u *)&mat(i + 4, j + 2),
                           *(__m128i *)&resultr3r);
          _mm_storeu_ps(&mat(i, j + 3), resultr4l);
          _mm_storel_epi64((__m128i_u *)&mat(i + 4, j + 3),
                           *(__m128i *)&resultr4r);
          _mm_storeu_ps(&mat(i, j + 4), resultr5l);
          _mm_storel_epi64((__m128i_u *)&mat(i + 4, j + 4),
                           *(__m128i *)&resultr5r);
          _mm_storeu_ps(&mat(i, j + 5), resultr6l);
          _mm_storel_epi64((__m128i_u *)&mat(i + 4, j + 5),
                           *(__m128i *)&resultr6r);
        }
        for (; j < YLENGTH - 1; ++j) {
          for (int ii = i; ii < i + 6; ++ii) {
            mat(ii, j) = (dy * dy * (buffer(ii + 1, j) + buffer(ii - 1, j)) +
                          dx * dx * (buffer(ii, j + 1) + buffer(ii, j - 1))) /
                         (2 * (dx * dx + dy * dy));
          }
        }
        for (; i + 6 >= XLENGTH - 1 - 8 && i < XLENGTH - 1; ++i) {
          for (int j = 1; j < YLENGTH - 1; ++j) {
            mat(i, j) = (dy * dy * (buffer(i + 1, j) + buffer(i - 1, j)) +
                         dx * dx * (buffer(i, j + 1) + buffer(i, j - 1))) /
                        (2 * (dx * dx + dy * dy));
          }
        }
      }
    } // end buffer -> mat

    refresh(a, mat, K);
  }
  std::clog << "离开calculate_sse()\n";
}


void calculate_cpu_buffer(txtRead &a, Eigen::MatrixXf &mat, int iters = 1,
                          double K = -10.) {
  double xmax = XMAX, xmin = XMIN;
  double ymax = YMAX, ymin = YMIN;
  double dx = (xmax - xmin) / XLENGTH;
  double dy = (ymax - ymin) / YLENGTH; // 空间步长
  auto buffer = mat;
  refresh(a, mat, K);
  refresh(a, buffer, K);
  for (int n = 0; n < iters; ++n) {
#pragma omp parallel for num_threads(8)
    for (auto i = 1; i < XLENGTH - 1; i++) {
      for (auto j = 1; j < YLENGTH - 1; j++) {
        buffer(i, j) = (dy * dy * (mat(i + 1, j) + mat(i - 1, j)) +
                        dx * dx * (mat(i, j + 1) + mat(i, j - 1))) /
                       (2 * (dx * dx + dy * dy));
      }
    }
    refresh(a, buffer, K);
#pragma omp parallel for num_threads(8)
    for (auto i = 1; i < XLENGTH - 1; i++) {
      for (auto j = 1; j < YLENGTH - 1; j++) {
        mat(i, j) = (dy * dy * (buffer(i + 1, j) + buffer(i - 1, j)) +
                     dx * dx * (buffer(i, j + 1) + buffer(i, j - 1))) /
                    (2 * (dx * dx + dy * dy));
      }
    }
    refresh(a, mat, K);
  }
}
void calculate_cpu(txtRead &a, Eigen::MatrixXf &mat, int iters = 1,
                   double K = -10.) {
  double xmax = XMAX, xmin = XMIN;
  double ymax = YMAX, ymin = YMIN;
  double dx = (xmax - xmin) / XLENGTH;
  double dy = (ymax - ymin) / YLENGTH; // 空间步长
  refresh(a, mat, K);
  for (int n = 0; n < iters; ++n) {
    double err = 0;
    // #pragma omp parallel for num_threads(12)
    for (auto i = 1; i < XLENGTH - 1; i++) {
      for (auto j = 1; j < YLENGTH - 1; j++) {
        float old = mat(i, j);
        mat(i, j) = (dy * dy * (mat(i + 1, j) + mat(i - 1, j)) +
                     dx * dx * (mat(i, j + 1) + mat(i, j - 1))) /
                    (2 * (dx * dx + dy * dy));
        err += std::abs(mat(i, j) - old);
      }
    }
    // 打印每次迭代都相对于上一次变化了多少，
    // 用这个i变化值除以l面积，y当作误差值，
    // 研究迭代次数与误差值之间的关系，并找到合适的阈值以停止迭代
    err /= XLENGTH * (float)YLENGTH;
    if (n % (12 * 10) == 0)
      std::cout << "\r                                     \r"
                << "iters: " << n << "\terr:" << err << std::flush;
    if (n == iters - 1)
      std::cout << "\n";
    refresh(a, mat, K);
  }
}
void calculate(txtRead &a, Eigen::MatrixXf &mat, int iters = 1,
               double K = -10.) {
  // 先遍历四周，然后遍历内部
  // 若矩阵空间为笛卡尔系，原点为矩阵左下，横轴为x，纵轴为y。
  // tip:下由于恒为1，不用给
  // 计算值
  calculate_sse(a, mat, iters, K);

  // 中间
  /*
  for(auto i=1;i<XLENGTH-1;i++){
      for(auto j=1;j<YLENGTH-1;j++){
          mat(i,j)=(dy*dy*(mat(i+1,j)+mat(i-1,j))+dx*dx*(mat(i,j+1)+mat(i,j-1)))/(2*(dx*dx+dy*dy));
      }
  }
  */
  /*
  for_ysh(1,XLENGTH-2,[&](int i){
      for_ysh(1,YLENGTH-2,[&](int j){
          mat(i,j)=(dy*dy*(mat(i+1,j)+mat(i-1,j))+dx*dx*(mat(i,j+1)+mat(i,j-1)))/(2*(dx*dx+dy*dy));
      });
  });*/
  /*
#pragma omp parallel for num_threads(8)
for(int i = 1; i < XLENGTH-1; ++i) {
  for(int j = 1; j < YLENGTH-1; ++j) {
      buffer(i, j) =
(dy*dy*(mat(i+1,j)+mat(i-1,j))+dx*dx*(mat(i,j+1)+mat(i,j-1)))/(2*(dx*dx+dy*dy));
  }
}
*/
}

void calculate_cap(txtRead &a, Eigen::MatrixXf &mat, DieManager &die,
                   double K = -10.) {
  double sum{0};
  double sum_cap(txtRead & a, Eigen::MatrixXf & mat, DieManager & die,
                 std::string name, double K);
  // 只有mental
  for (auto i : a.mental) {
    if (i.get_name() == "c2") {
      sum = sum_cap(a, mat, die, "c2", K);
      std::cout << "c2=" << sum << "\n";
    }
    //std::cout << "\n\n\n\n\n";
    if (i.get_name() == "rEnv") {
      sum = sum_cap(a, mat, die, "rEnv", K);
      std::cout << "rEnv=" << sum << "\n";
    }
  }
}

double sum_cap(txtRead &a, Eigen::MatrixXf &mat, DieManager &die,
               std::string name, double K = -10.) {
  double xmax = XMAX, xmin = XMIN;
  double ymax = YMAX, ymin = YMIN;
  double dx = (xmax - xmin) / XLENGTH;
  double dy = (ymax - ymin) / YLENGTH; // 空间步长
  float sum{0};
  for (auto i : a.mental) {
    if (i.get_name() != name) {
      continue;
    }
    // 左右
    double max{0};
    for (auto j : i)
      if (j.y() > max)
        max = j.y();
    for (int y_index = 0, max_y_index = (int)std::ceil((max - i[0].y()) / dy);
         y_index < max_y_index; ++y_index) {
      float y = i[0].y() + y_index * dy; // 当前循环y实际值
      float xl = i[0].x() + y_index * dy / K - dx;
      float xr = i[1].x() - y_index * dy / K + dx;
      int xl_index = (xl - xmin) / dx;
      int xr_index = (xr - xmin) / dx;
      int indeed_y_index = (y - ymin) / dy;
      { // 左
        auto result = die.search_die({xl, y});
        float unwrap_result;
        if (result) {
          unwrap_result = *result;
        } else {
          std::cerr << "查找到非法点1\n";
          std::cout << xl << " " << y << "\n";
          exit(0);
        }
        sum = sum + mat(xl_index, indeed_y_index) * unwrap_result * dy / dx;
      }
      { // 右，加一判断
        auto result = die.search_die({xr, y});
        float unwrap_result;
        if (result) {
          unwrap_result = *result;
        } else {
          std::cerr << "查找到非法点2\n";
          std::cout << xr << " " << y << "\n";
          exit(0);
        }
        if (mat(xr_index - 1, indeed_y_index) != 0) {
          sum =
              sum + mat(xr_index - 1, indeed_y_index) * unwrap_result * dy / dx;
        } else {
          sum = sum + mat(xr_index, indeed_y_index) * unwrap_result * dy / dx;
        }
      }
    }
    // 上&下
    // 下
    {
      for (int x_index = 0,
               max_x_index = (int)std::ceil((i[1].x() - i[0].x()) / dx);
           x_index < max_x_index; ++x_index) {
        double y = i[0].y() - dy;
        int y_index = (y - ymin) / dy;
        float x = x_index * dx + i[0].x();
        int xd_index = (x - xmin) / dx;
        {
          auto result = die.search_die({x, y});
          float unwrap_result;
          if (result) {
            unwrap_result = *result;
          } else {
            std::cerr << "查找到非法点3\n";
            std::cout << x << " " << y << "\n";
            exit(0);
          }
          sum = sum + mat(xd_index, y_index) * unwrap_result * dx / dy;
        }
      }
    }
    // 上
    {
      double xmax = -10;
      for (auto j : i)
        if (j.x() > xmax)
          xmax = j.x();
      double ymax = 0;
      for (auto j : i)
        if (j.y() > ymax)
          ymax = j.y();
      for (int x_index = 0,
               max_x_index = (int)std::ceil(
                   (i[1].x() - i[0].x() + 2 * (xmax - i[1].x())) / dx);
           x_index < max_x_index; ++x_index) {
        double y = ymax + dy;
        int y_index = (y - ymin) / dy;
        float x = x_index * dx + i[0].x() - (xmax - i[1].x());
        int xd_index = (x - xmin) / dx;
        {
          auto result = die.search_die({x, y});
          float unwrap_result;
          if (result) {
            unwrap_result = *result;
          } else {
            std::cerr << "查找到非法点3\n";
            std::cout << x << " " << y << "\n";
            exit(0);
          }
          if (mat(xd_index, y_index - 1) != 0) {
            sum = sum + mat(xd_index, y_index - 1) * unwrap_result * dx / dy;
          } else {
            sum = sum + mat(xd_index, y_index) * unwrap_result * dx / dy;
          }
        }
      }
    }
  }
  return sum;
}
