#include <cstddef>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <txtRead/txtRead.hpp>
#include <vector>
#include <CLI/CLI.hpp>

extern void copredict1(const std::vector<double>&, std::ostream& = std::cout);
extern void copredict2(const std::vector<double>&, std::ostream& = std::cout);
extern void copredict3(const std::vector<double>&, std::ostream& = std::cout);
extern void copredict4(const std::vector<double>&, std::ostream& = std::cout);
extern void copredict5(const std::vector<double>&, std::ostream& = std::cout);
extern void copredict6(const std::vector<double>&, std::ostream& = std::cout);
extern void copredict7(const std::vector<double>&, std::ostream& = std::cout);

int main(int argc, char *argv[]) {
  CLI::App app{"solver"};

  std::string outputFilename{};
  std::string inputFilename{};

  app.add_option("-o,--output", outputFilename, "Output filename");
  app.add_option("-i,--input", inputFilename, "Input filename");

  CLI11_PARSE(app, argc, argv);

  if (outputFilename.empty() || inputFilename.empty()) {
    std::cerr << "Error: Both -o and -i options are required." << std::endl;
    return 1;
  }

  txtRead a;
  double w1 = -1, w2 = -2, edge = -3, edge_space = -4;
  a.Open(inputFilename);
  a.Read();
  std::vector<double> input_data{};
  // std::function<void()> funs[] = {f1, f2};
  // int n = 0;
  // funs[n]();
  cell master{a.get_master()};
  CellVec mental{a.get_mental()};
  CellVecAdd die{a.get_die()};
  cell corner{a.get_corner()};

  // 判断是否为patern1
  double max = 0;
  std::vector<double> max_value{};

  for (auto i : master) {
    if (i.y() > max) {
      max = i.y();
    }
  }
  max_value.push_back(max);

  for (auto i : mental) {
    max = 0;
    for (auto j : i) {
      if (j.y() > max) {
        max = j.y();
      }
    }

    for (auto j : max_value) {
      if (j == max) {
        break;
      }
      if (j == max_value.back()) {
        max_value.push_back(max);
      }
    }
  }
  // std::cout << max_value.size() << "\n";
  // 操作结束，此时max_value.size()为2则表示pattern1，其他为3,然后读取corner判断pattern2与pattern3
  if(std::fabs(master.data()[1].x()-master.data()[2].x())>1e-6||std::fabs(mental[2][1].x()-mental[2][2].x())>1e-6||std::fabs(mental[1][1].x()-mental[1][2].x())>1e-6){
      //pattern123
      if(max_value.size()==2){
          input_data.push_back(max_value[0]);
      }else{
        std::sort(std::begin(max_value), std::end(max_value));
        input_data.push_back(max_value[1]);
        input_data.push_back(max_value[2]);
      }

      if (max_value.size() == 2) {
        input_data.push_back(1.);
        goto pattern1;
      } else {
        // 判断2 or 3
        max = 0;
        for (auto i : corner) {
          if (i.y() > max) {
            max = i.y();
          }
        }
        if (max >= 6) {
          input_data.push_back(3.);
          goto pattern3;
        } else {
          input_data.push_back(2.);
          goto pattern2;
        }
      }
  }else{
      //pattern4567
      //p4=diag p5=plate2l p6=plate4l p7=stack3l
      //先存入高度
      if(max_value.size()==2){
          input_data.push_back(max_value[0]);
      }else{
        std::sort(std::begin(max_value), std::end(max_value));
        input_data.push_back(max_value[1]);
        input_data.push_back(max_value[2]);
      }

    if(max_value.size()==2){
        input_data.push_back(5.);
        goto pattern5;
    }else{
        max = 0;
        for (auto i : corner) {
          if (i.y() > max) {
            max = i.y();
          }
        }
        if (max <= 6) {
          input_data.pop_back();
          input_data.push_back(6.);
          goto pattern6;
        } else {
            for(auto &i:mental){
                if(i.get_name()!="botleft" && i.get_name()!="botright"){
                    for(auto &j:i){
                        if(std::fabs(j.x()-0)<1e-6){
                            input_data.push_back(7.);
                            goto pattern7;
                        }
                    }
                }
            }
            input_data.push_back(4.);
            goto pattern4;
        }

    }
  }
  // pattern1
pattern1:
  //std::cout << "this is pattern1\n";
  if (mental.get_number() == 2) {
    edge = 4.5;
    edge_space = 4.5;

    max = 0;
    for (auto i : master) {
      if (i.x() > max) {
        max = i.x();
      }
    }
    w1 = ((master[1].x() - master[0].x())) / 2 + max;
    max = 0;
    for (auto i : corner) {
      if (i.x() > max) {
        max = i.x();
      }
    }
    w2 = (max - 0.5 - 9 - w1 / 2) / 2;
    goto output; // 输出w1,w2,edge,edge_space
  } else if (mental.get_number() == 3) {
    // 判断有两个导体的情况（不考虑下方导体），有master-mental[2](edge==4.5),master-mental[3](edge_space==4.5)两种情况
    int count = 0;
    for (auto i : mental) {
      if (i.get_name() == "c2") {
        count = 1;
      }
      if (i.get_name() == "lEnv") {
        count = 2;
      }
    }
    if (count == 1) {
      // edge_space==4.5
      edge_space = 4.5;
      max = 0;
      for (auto i : master) {
        if (i.x() > max) {
          max = i.x();
        }
      }
      w1 = ((master[1].x() - master[0].x())) / 2 + max;

      max = 0;
      for (auto i : mental[2]) {
        if (i.x() > max) {
          max = i.x();
        }
      }
      w2 = (mental[2][1].x() - mental[2][0].x()) / 2 + max -
           (mental[2][1].x() + mental[2][0].x()) / 2;
      edge = (mental[2][0].x() - master[0].x()) - w1;
      goto output;
    }

    if (count == 2) {
      // edge==4.5
      edge = 4.5;

      max = 0;
      for (auto i : master) {
        if (i.x() > max) {
          max = i.x();
        }
      }
      w1 = ((master[1].x() - master[0].x())) / 2 + max;
      edge_space = (master[0].x() - mental[2][0].x()) - w1;

      max = 0;
      for (auto i : corner) {
        if (i.x() > max) {
          max = i.x();
        }
      }
      w2 = (max - 0.5 - w1 / 2 - edge - edge_space) / 2;

      goto output;
    }
  } else {
    // 四个导体情况
    max = 0;
    for (auto i : master) {
      if (i.x() > max) {
        max = i.x();
      }
    }
    w1 = (master[1].x() - master[0].x()) / 2 + max;

    max = 0;
    for (auto i : mental[3]) {
      if (i.x() > max) {
        max = i.x();
      }
    }
    w2 = (mental[3][1].x() - mental[3][0].x()) / 2 + max -
         (mental[3][1].x() + mental[3][0].x()) / 2;

    edge = mental[3][0].x() - master[0].x() - w1;
    edge_space = master[0].x() - mental[2][0].x() - w1;
    goto output;
  }

// pattern2
pattern2:
  //std::cout << "this is pattern2\n";
  // 共分为4，5，7，9个导体的情况，4为均4.5的12文件，5为e=4.5的12文件u与e_s=4.5的1文件，
  // 7为e_s=4.5的2文件和均无4.5的1文件，9为均无4.5的2文件
  // 4
  if (mental.get_number() == 4) {
    max = 0;
    for (auto i : master) {
      if (i.y() > max) {
        max = i.y();
      }
    }
    if (max < 0.1) {
      // 文件2
      edge = 4.5;
      edge_space = 4.5;
      max = 0;
      for (auto i : mental[1]) {
        if (i.x() > max) {
          max = i.x();
        }
      }
      w1 = (mental[1][1].x() - mental[1][0].x()) / 2 + max;
      w2 = (master[1].x() - edge) * 2 - w1;
      goto output;
    } else {
      // 文件1
      edge = 4.5;
      edge_space = 4.5;
      max = 0;
      for (auto i : master)
        if (i.x() > max)
          max = i.x();
      w1 = master[1].x() + max;
      max = 0;
      for (auto i : corner)
        if (i.x() > max)
          max = i.x();

      w2 = (max - 9.5 - w1 / 2) / 2;
      goto output;
    }
  }
  // 5
  if (mental.get_number() == 5) {
    if (mental[2].get_name() == "c2") {
      // edge_space=4.5 && file 1
      edge_space = 4.5;
      max = 0;
      for (auto i : master)
        if (i.x() > max)
          max = i.x();
      w1 = master[1].x() + max;
      max = 0;
      for (auto i : mental[2])
        if (i.x() > max)
          max = i.x();
      w2 = ((mental[2][1].x() - mental[2][0].x()) / 2 + max -
            (mental[2][1].x() + mental[2][0].x()) / 2);
      edge = mental[2][0].x() - master[0].x() - w1;
      goto output;
    } else {
      // edge=4.5 file 1/2
      max = 0;
      for (auto i : master)
        if (i.y() > max)
          max = i.y();
      if (max < 0.1) {
        // file2
        edge = 4.5;
        max = 0;
        for (auto i : mental[2])
          if (i.x() > max)
            max = i.x();
        w1 = (mental[2][1].x() - mental[2][0].x()) / 2 + max;
        edge_space = mental[2][0].x() - mental[1][0].x() - w1;
        max = 0;
        for (auto i : corner)
          if (i.x() > max)
            max = i.x();
        w2 = (max - 5 - w1 / 2 - edge_space) / 2;
        goto output;
      } else {
        // file1
        edge = 4.5;
        max = 0;
        for (auto i : master)
          if (i.x() > max)
            max = i.x();
        w1 = (master[1].x() - master[0].x()) / 2 + max;
        edge_space = master[0].x() - mental[2][0].x() - w1;
        max = 0;
        for (auto i : corner)
          if (i.x() > max)
            max = i.x();
        w2 = (max - 5 - w1 / 2 - edge_space) / 2;
        goto output;
      }
    }
  }
  // 7
  if (mental.get_number() == 7) {
    // 共有edge_space=4.5 file2 与均无4.5 file1两种情况
    max = 0;
    for (auto i : master)
      if (i.y() > max)
        max = i.y();
    if (max < 0.1) {
      // file2
      edge_space = 4.5;
      max = 0;
      for (auto i : mental[2])
        if (i.x() > max)
          max = i.x();
      w1 = (mental[2][1].x() - mental[2][0].x()) / 2 + max;
      edge = mental[3][0].x() - mental[2][0].x() - w1;
      max = 0;
      for (auto i : corner)
        if (i.x() > max)
          max = i.x();
      w2 = (max - 5 - edge - w1 / 2) / 2;
      goto output;
    } else {
      // file1
      max = 0;
      for (auto i : master)
        if (i.x() > max)
          max = i.x();
      w1 = (master[1].x() - master[0].x()) / 2 + max;
      max = 0;
      for (auto i : mental[3])
        if (i.x() > max)
          max = i.x();
      w2 = (mental[3][1].x() - mental[3][0].x()) / 2 + max -
           (mental[3][1].x() + mental[3][0].x()) / 2;
      edge = mental[3][0].x() - master[0].x() - w1;
      edge_space = master[0].x() - mental[2][0].x() - w1;
      goto output;
    }
  }

  // 9
  if (mental.get_number() == 9) {
    // 最后一种情况，均无4.5，file2
    max = 0;
    for (auto i : mental[3])
      if (i.x() > max)
        max = i.x();
    w1 = (mental[3][1].x() - mental[3][0].x()) / 2 + max;
    max = 0;
    for (auto i : mental[4])
      if (i.x() > max)
        max = i.x();
    w2 = (mental[4][1].x() - mental[4][0].x()) / 2 + max -
         (mental[4][1].x() + mental[4][0].x()) / 2;
    edge = mental[4][0].x() - mental[3][0].x() - w1;
    edge_space = mental[3][0].x() - mental[2][0].x() - w1;
    goto output;
  }

// pattern3
pattern3:
  //std::cout << "this is pattern3\n";
  // mental共有4，6，10三种，4对应均为4.5情况，6对应edge为4.5或edge_space为4.5情况，10对应均不为4.5情况
  // 4
  if (mental.get_number() == 4) {
    // 包括file12
    if (master.get_name() == "c1") {
      // file1
      edge = 4.5;
      edge_space = 4.5;
      max = 0;
      for (auto i : master)
        if (i.x() > max)
          max = i.x();
      w1 = (master[1].x() - master[0].x()) / 2 + max;

      max = 0;
      for (auto i : corner)
        if (i.x() > max)
          max = i.x();
      w2 = (max - 9.5 - w1 / 2) / 2;
      goto output;
    } else {
      // file2
      edge = 4.5;
      edge_space = 4.5;
      max = 0;
      for (auto i : mental[3])
        if (i.x() > max)
          max = i.x();
      w1 = mental[3][1].x() + max;
      max = 0;
      for (auto i : corner)
        if (i.x() > max)
          max = i.x();
      w2 = (max - 9.5 - w1 / 2) / 2;
      goto output;
    }
  }
  // 6
  if (mental.get_number() == 6) {
    // 共有四个，分别为edge=4.5,file12ledge_space=4.5,file12
    if (master.get_name() == "c1") {
      // 为两个file1
      max = 0;
      for (auto i : master)
        if (i.x() > max)
          max = i.x();
      w1 = master[1].x() + max;

      if (mental[2].get_name() == "c3") {
        // edge_space=4.5
        edge_space = 4.5;
        max = 0;
        for (auto i : mental[2])
          if (i.x() > max)
            max = i.x();
        w2 = (mental[2][1].x() - mental[2][0].x()) / 2 + max -
             (mental[2][1].x() + mental[2][0].x()) / 2;
        edge = mental[2][0].x() - master[0].x() - w1;
        goto output;
      } else {
        // edge=4.5
        edge = 4.5;
        edge_space = master[0].x() - mental[2][0].x() - w1;
        max = 0;
        for (auto i : corner)
          if (i.x() > max)
            max = i.x();
        w2 = (max - 0.5 - edge - edge_space - w1 / 2) / 2;
        goto output;
      }
    } else {
      // 为两个file2
      if (mental[4].get_name() == "c3") {
        // edge_space=4.5
        edge_space = 4.5;
        max = 0;
        for (auto i : mental[3])
          if (i.x() > max)
            max = i.x();
        w1 = mental[3][1].x() + max;
        max = 0;
        for (auto i : mental[4])
          if (i.x() > max)
            max = i.x();
        w2 = (mental[4][1].x() - mental[4][0].x()) / 2 + max -
             (mental[4][1].x() + mental[4][0].x()) / 2;
        edge = mental[4][0].x() - mental[2][0].x() - w1;

        goto output;
      } else {
        // edge=4.5
        edge = 4.5;
        max = 0;
        for (auto i : mental[4])
          if (i.x() > max)
            max = i.x();
        w1 = mental[4][1].x() + max;
        edge_space = mental[3][0].x() - mental[2][0].x() - w1;
        max = 0;
        for (auto i : corner)
          if (i.x() > max)
            max = i.x();
        w2 = (max - 0.5 - edge - edge_space - w1 / 2) / 2;
        goto output;
      }
    }
  }
  // 10
  if (mental.get_number() == 10) {
    // 包括file12
    if (master.get_name() == "c1") {
      // file1
      max = 0;
      for (auto i : master)
        if (i.x() > max)
          max = i.x();
      w1 = master[1].x() + max;
      max = 0;
      for (auto i : mental[3])
        if (i.x() > max)
          max = i.x();
      w2 = (mental[3][1].x() - mental[3][0].x()) / 2 + max -
           (mental[3][1].x() + mental[3][0].x()) / 2;
      edge_space = master[0].x() - mental[2][0].x() - w1;
      edge = mental[3][0].x() - master[0].x() - w1;
      goto output;
    } else {
      // file2
      max = 0;
      for (auto i : mental[4])
        if (i.x() > max)
          max = i.x();
      w1 = mental[4][1].x() + max;
      max = 0;
      for (auto i : mental[5])
        if (i.x() > max)
          max = i.x();
      w2 = (mental[5][1].x() - mental[5][0].x()) / 2 + max -
           (mental[5][1].x() + mental[5][0].x()) / 2;
      edge = mental[5][0].x() - mental[3][0].x() - w1;
      edge_space = mental[3][0].x() - mental[2][0].x() - w1;

      goto output;
    }
  }
pattern4:
  //std::cout<<"This is pattern4 \n";
  //变量名->实际意义  w1->w1  w2->edgespace  edge->s1  edge_space->s2
  w1=master[1].x()-master[0].x();
  w2=master[0].x()-mental[2][1].x();
  edge=mental[3][0].x()-master[1].x();
  edge_space=mental[6][0].x()-mental[5][1].x();
  goto output;
pattern5:
  //std::cout<<"This is pattern5 \n";
  //plate2l  w1->w1  w2->s1 edge->s2  edge_space->1  
  w1=master[1].x()-master[0].x();
  w2=master[0].x()-mental[2][1].x();
  edge=mental[3][0].x()-master[1].x();
  edge_space=1;
  goto output;
pattern6:
  //std::cout<<"This is pattern6 \n";
  //plate3l  w1->w1  w2->s1  edge->s2  edge_space->1
  w1=master[1].x()-master[0].x();
  w2=master[0].x()-mental[2][1].x();
  edge=mental[3][0].x()-master[1].x();
  edge_space=1;
  goto output;
pattern7:
  //std::cout<<"This is pattern7 \n";
  //stack3l  w1->w1  w2->s1  edge->s2  edge_space->w3
  w1=master[1].x()-master[0].x();
  w2=master[0].x()-mental[2][1].x();
  edge=mental[3][0].x()-master[1].x();
  edge_space=mental[6][1].x()-mental[6][0].x();
  goto output;


  // 输出
output:
  double pattern_switch = input_data.back();
  input_data.pop_back();
  input_data.push_back(w1);
  input_data.push_back(w2);
  input_data.push_back(edge);
  input_data.push_back(edge_space);
  std::ofstream ofs{outputFilename};
  if(std::abs(pattern_switch - 1.) < 1e-6)
    copredict1(input_data, ofs);
  else if(std::abs(pattern_switch - 2.) < 1e-6)
    copredict2(input_data, ofs);
  else if(std::abs(pattern_switch - 3.) < 1e-6)
    copredict3(input_data, ofs);
  else if(std::abs(pattern_switch - 4.) < 1e-6)
    copredict4(input_data, ofs);
  else if(std::abs(pattern_switch - 5.) < 1e-6)
    copredict5(input_data, ofs);
  else if(std::abs(pattern_switch - 6.) < 1e-6)
      copredict6(input_data, ofs);
  else if(std::abs(pattern_switch - 7.) < 1e-6)
    copredict7(input_data, ofs);

}
