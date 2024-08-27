#include "./include/fem/full_matrix.h"
#include "./include/fem/full_vector.h"
#include "./include/fem/base_function.h"
#include <iostream>

fem::FullVector conjugate_gradient(const fem::FullMatrix &A,
                                   const fem::FullVector &b, int max_iter,
                                   double tolerance) {
  std::cerr << "设置初始值\n";
  fem::FullVector x(b.rows());
  fem::FullVector r = A * x - b;
  fem::FullVector p = r;

  double rsold = r * r;

  std::cerr << "进入循环\n";
  for (int iter = 0; iter < max_iter; ++iter) {
    fem::FullVector Ap = A * p;
    double alpha = rsold / ((p * Ap) + 1e-10);

    x = x + p * alpha;
    r = r - Ap * alpha;

    double rsnew = r * r;

    if (r.norm() < tolerance) {
      std::cout << "Convergence achieved at iteration " << iter << "\n";
      break;
    }
    if (iter == max_iter - 1) {
      std::cout << "max iter\n";
    }

    p = r + p * (rsnew / rsold);
    rsold = rsnew;
  }
  return x;
}

fem::FullVector jacobi_iteration(const fem::FullMatrix &A,
                                 const fem::FullVector &b, double tol,
                                 int max_iter) {
  int n = b.rows();
  fem::FullVector x(b.rows());

  for (int k = 0; k < max_iter; ++k) {
    fem::FullVector x_old = x;
    for (int i = 0; i < n; ++i) {
      double sigma = 0;
      for (int j = 0; j < n; ++j) {
        if (j != i) {
          sigma += A(i, j) * x(j, 0);
        }
      }
      x(i, 0) = (b(i, 0) - sigma) / A(i, i);
    }
    if (fem::FullVector(x - x_old).norm() < tol) {
      std::cout << "Convergence achieved at iteration " << k << "\n";
      return x;
    }
  }
  std::cout << "Warning: The iteration did not converge after " << max_iter
            << " iterations.\n";
  return x;
}
int main() {
    std::cout << fem::AttenuationFunction1D<double>(0.6) << "\n";
    return 0;
  fem::FullMatrix a(2, 2);
  fem::FullVector b(2);
  a(0, 0) = 100;
  a(0, 1) = 10;
  a(1, 0) = 1000;
  a(1, 1) = 10000;

  b(0, 0) = 100*0.618 + 10 * 3.1415926;
  b(1, 0) = 1000 * 0.618 + 10000* 3.1415926;
  /*
  a(0, 0) = 10;
  a(0, 1) = -1;
  a(0, 2) = 2;
  a(0, 3) = 0;
  a(1, 0) = -1;
  a(1, 1) = 11;
  a(1, 2) = -1;
  a(1, 3) = 3;
  a(2, 0) = 2;
  a(2, 1) = -1;
  a(2, 2) = 10;
  a(2, 3) = -1;
  a(3, 0) = 0;
  a(3, 1) = 3;
  a(3, 2) = -1;
  a(3, 3) = 8;

  b(0, 0) = 6;
  b(1, 0) = 25;
  b(2, 0) = -11;
  b(3, 0) = 15;
  */
  try {
    std::cerr << "开始迭代\n";
    // auto c = conjugate_gradient(a, b, 1000, 1e-6);
    auto c = jacobi_iteration(a, b, 1e-6, 1000);

    std::cout << "A:\n";
    a.print(std::cout);
    std::cout << "b:\n";
    b.print(std::cout);
    std::cout << "x:\n";
    c.print(std::cout);
  } catch (const char *e) {
    std::cerr << e;
  }
}
