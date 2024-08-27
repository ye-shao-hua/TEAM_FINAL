#pragma once
template <class FUNC> void for_ysh(int begin, int end, FUNC func) {
#pragma omp parallel for num_threads(8)
  for (auto i = begin; i <= end; i++) {
    func(i);
  }
}
