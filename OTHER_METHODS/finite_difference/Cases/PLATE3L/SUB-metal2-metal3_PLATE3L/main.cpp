#include <fstream>
#include <string>

int main() {
  std::ifstream ifs{"SUB-metal2-metal3_PLATE3L.tbl.text"};
  std::ofstream ofs{"p2m23.text"};
  while (1) {
    std::string buffer;
    std::getline(ifs, buffer);
    if (ifs.eof()) {
      break;
    }
    ofs << 0.4665 + 0.187 << "\t " << 0.4665 + 2 * 0.187 << "\t " << buffer
        << "\n";
  }
}
