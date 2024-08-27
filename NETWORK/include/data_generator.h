#include <fstream>
#include <sstream>
#include <string_view>
#include <tiny_dnn/tiny_dnn.h>
#include <utility>
#include <vector>

namespace data_generator {
using namespace tiny_dnn;
using namespace std;

using input_data = vector<vec_t>;
using output_data = vector<vec_t>;

pair<input_data, output_data> read_text(string_view filename) {
  ifstream ifs{filename.data()};
  string buffer;
  getline(ifs, buffer);
  int number_of_line = 0;
  int split_index = 0;

  for (istringstream iss{buffer}; !iss.eof(); ++number_of_line) {
    iss >> buffer;
    if (buffer == "|")
      split_index = number_of_line;
  }
  std::cout << number_of_line << ' ' << split_index << '\n';

  vector<vec_t> input_data, output_data;

  do {
    getline(ifs, buffer);
    istringstream iss{buffer};
    vec_t line_input_data;
    for (int i = 0; i < split_index; ++i) {
      double number = 0;
      iss >> number;
      line_input_data.push_back(number);
    }
    if (line_input_data[0] == 0)
      break;
    input_data.push_back(line_input_data);

    string ignore;
    iss >> ignore;

    vec_t line_output_data;
    for (int i = split_index + 1; i < number_of_line; ++i) {
      double number;
      iss >> number;
      line_output_data.push_back(number);
    }
    output_data.push_back(line_output_data);
  } while (!ifs.eof());
  return {input_data, output_data};
}

} // namespace data_generator
