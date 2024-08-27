#include <vector>
#include "./full_matrix.h"
#include "./point.h"
namespace fem{
    template<int dim,typename NumberType>
    class MappingArray{
        public:
            MappingArray(Point<dim,NumberType> *point);
        private:
            std::vector<FullMatrix> _data;
    };
}
