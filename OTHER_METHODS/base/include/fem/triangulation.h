#include "./cell.h"
#include <vector>
#include "./hyper_cube.h"
namespace fem{
    template<int dim ,typename NumberType>
    class Triangulation{
        public:
        Triangulation(HyperCube<dim,NumberType> &hy);
        void refine(int id);
        private:
            std::vector<Cell<dim,NumberType>> _data;
    };

}
