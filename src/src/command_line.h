// Copyright (c) 2015, The Regents of the University of California (Regents)
// See LICENSE.txt for license details

#ifndef COMMAND_LINE_H_
#define COMMAND_LINE_H_

#include <getopt.h>
#include <libgen.h>

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>


/*
GAP Benchmark Suite
Class:  CLBase
Author: Scott Beamer

Handles command line argument parsing
 - Through inheritance, can add more options to object
 - For example, most kernels will use CLApp
*/


class CLBase {
 protected:
  int argc_;
  char** argv_;
  std::string name_;
  std::string get_args_ = "f:g:hk:su:q:";
  std::vector<std::string> help_strings_;

  int scale_ = -1;
  int degree_ = 16;
  std::string filename_ = "";
  bool symmetrize_ = false;
  bool uniform_ = false;

  void AddHelpLine(char opt, std::string opt_arg, std::string text,
                   std::string def = "") {
    const int kBufLen = 100;
    char buf[kBufLen];
    if (opt_arg != "")
      opt_arg = "<" + opt_arg + ">";
    if (def != "")
      def = "[" + def + "]";
    snprintf(buf, kBufLen, " -%c %-9s: %-54s%10s", opt, opt_arg.c_str(),
            text.c_str(), def.c_str());
    help_strings_.push_back(buf);
  }

 public:
  CLBase(int argc, char** argv, std::string name = "") :
         argc_(argc), argv_(argv), name_(name) {
    AddHelpLine('h', "", "print this help message");
    AddHelpLine('f', "file", "load graph from file");
    AddHelpLine('s', "", "symmetrize input edge list", "false");
    AddHelpLine('g', "scale", "generate 2^scale kronecker graph");
    AddHelpLine('u', "scale", "generate 2^scale uniform-random graph");
    AddHelpLine('k', "degree", "average degree for synthetic graph",
                std::to_string(degree_));
    AddHelpLine('q', "threads", "set number of threads");
  }

  bool ParseArgs() {
    signed char c_opt;
    extern char *optarg;          // from and for getopt
    while ((c_opt = getopt(argc_, argv_, get_args_.c_str())) != -1) {
      HandleArg(c_opt, optarg);
    }
    if ((filename_ == "") && (scale_ == -1)) {
      std::cout << "No graph input specified. (Use -h for help)" << std::endl;
      return false;
    }
    if (scale_ != -1)
      symmetrize_ = true;
    return true;
  }

  std::string getGraphBasename() {
      if (filename_ != "") return std::string(basename(((char*)(filename_.c_str()))));
      if (uniform_) return "uniform";
      if (scale_ != -1) return "kronecker";
      return "unknown";
  }

  int getThreadNum() {
      auto value = getenv("OMP_NUM_THREADS");
      if (value == NULL) return 1;
      return std::stoi(value);
  }

  void virtual HandleArg(signed char opt, char* opt_arg) {
    switch (opt) {
      case 'f': filename_ = std::string(opt_arg);           break;
      case 'g': scale_ = atoi(opt_arg);                     break;
      case 'h': PrintUsage();                               break;
      case 'k': degree_ = atoi(opt_arg);                    break;
      case 's': symmetrize_ = true;                         break;
      case 'u': uniform_ = true; scale_ = atoi(opt_arg);    break;
      case 'q': setenv("OMP_NUM_THREADS", opt_arg, 1);      break;
    }
  }

  void PrintUsage() {
    std::cout << name_ << std::endl;
    // std::sort(help_strings_.begin(), help_strings_.end());
    for (std::string h : help_strings_)
      std::cout << h << std::endl;
    std::exit(0);
  }

  int scale() const { return scale_; }
  int degree() const { return degree_; }
  std::string filename() const { return filename_; }
  bool symmetrize() const { return symmetrize_; }
  bool uniform() const { return uniform_; }
};



class CLApp : public CLBase {
  bool do_analysis_ = false;
  int num_trials_ = 1;
  int64_t start_vertex_ = -1;
  bool do_verify_ = false;
  float tau_ = 0.5;
  int m_ = 0;
  int k_ = 0;
  float p_ = 1.0;

 public:
 void set_verify() {do_verify_ = true;}
 
  CLApp(int argc, char** argv, std::string name) : CLBase(argc, argv, name) {
    get_args_ += "an:r:vy:m:b:p:";
    AddHelpLine('a', "", "output analysis of last run", "false");
    AddHelpLine('n', "n", "perform n trials", std::to_string(num_trials_));
    AddHelpLine('r', "node", "start from node r", "rand");
    AddHelpLine('v', "", "verify the output of each run", "false");
    AddHelpLine('y', "tau", "clustering measure", std::to_string(5));
    AddHelpLine('m', "m", "byte-array size for BF (default = 0 -> automatic)", std::to_string(0));
    AddHelpLine('b', "b", "hasing funtion k for BF (default = 0 -> automatic)", std::to_string(0));

  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'a': do_analysis_ = true;                    break;
      case 'n': num_trials_ = atoi(opt_arg);            break;
      case 'r': start_vertex_ = atol(opt_arg);          break;
      case 'v': do_verify_ = true;                      break;
      case 'y': tau_ = atof(opt_arg);                   break;
      case 'm': m_ = atoi(opt_arg);                   break;
      case 'b': k_ = atoi(opt_arg);                   break;
      case 'p': p_ = atof(opt_arg);                   break;
      default: CLBase::HandleArg(opt, opt_arg);
    }
  }

  bool do_analysis() const { return do_analysis_; }
  int num_trials() const { return num_trials_; }
  int64_t start_vertex() const { return start_vertex_; }
  bool do_verify() const { return do_verify_; }
  float tau() const { return tau_; }
  int m() const  {return m_;}
  int k() const {return k_;}
  float p() const {return p_;}
};



class CLIterApp : public CLApp {
  int num_iters_;

 public:
  CLIterApp(int argc, char** argv, std::string name, int num_iters) :
    CLApp(argc, argv, name), num_iters_(num_iters) {
    get_args_ += "i:";
    AddHelpLine('i', "i", "perform i iterations", std::to_string(num_iters_));
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'i': num_iters_ = atoi(opt_arg);            break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  int num_iters() const { return num_iters_; }
};

class CLSIMDApp: public CLApp {
  float treshold_ = 3;
  FILE *fp_tr = nullptr;

 public:
  CLSIMDApp(int argc, char** argv, std::string name) :
    CLApp(argc, argv, name){
    get_args_ += "t:";
    AddHelpLine('t', "t", "use simd treshold t", std::to_string(treshold_));
    get_args_ += "x:";
    AddHelpLine('x', "x", "file to load parameters (thresholds, representation etc)");
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 't': 
        treshold_ = atof(opt_arg);            
        break;
      case 'x': 
        fp_tr = fopen(opt_arg, "r");            
        if (fp_tr == nullptr) {
          printf("Could not open file %s\n", opt_arg);
          exit(-1);
        }
        break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  float treshold() const { return treshold_; }
  FILE* get_parameters_file() const { return fp_tr; }
};

class CLBlockApp: public CLApp {
  float k_scale_ = 0.1;
  int block_size_ = 1000;

 public:
  CLBlockApp(int argc, char** argv, std::string name) :
    CLApp(argc, argv, name){
    get_args_ += "z:c:";
    AddHelpLine('z', "z", "set blocksize", std::to_string(block_size_));
    AddHelpLine('c', "c", "fraction of blockisze used to sum up the block c", std::to_string(k_scale_));
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'c': k_scale_ = atof(opt_arg);            break;
      case 'z': block_size_ = atoi(opt_arg);            break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  float k_scale() const { return k_scale_; }
  float block_size() const { return block_size_; }
  
};

class CLBloomApp: public CLApp {
  float accuracy_ = 0.1;

 public:
  CLBloomApp(int argc, char** argv, std::string name) :
    CLApp(argc, argv, name){
    get_args_ += "c:";
    AddHelpLine('c', "c", "use bf accuracy c", std::to_string(accuracy_));
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'c': accuracy_ = atof(opt_arg);            break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  float accuracy() const { return accuracy_; }
};

class CLPageRank : public CLApp {
  int max_iters_;
  double tolerance_;

 public:
  CLPageRank(int argc, char** argv, std::string name, double tolerance,
             int max_iters) :
    CLApp(argc, argv, name), max_iters_(max_iters), tolerance_(tolerance) {
    get_args_ += "i:t:";
    AddHelpLine('i', "i", "perform at most i iterations",
                std::to_string(max_iters_));
    AddHelpLine('t', "t", "use tolerance t", std::to_string(tolerance_));
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'i': max_iters_ = atoi(opt_arg);            break;
      case 't': tolerance_ = std::stod(opt_arg);            break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  int max_iters() const { return max_iters_; }
  double tolerance() const { return tolerance_; }
};



template<typename WeightT_>
class CLDelta : public CLApp {
  WeightT_ delta_ = 1;

 public:
  CLDelta(int argc, char** argv, std::string name) : CLApp(argc, argv, name) {
    get_args_ += "d:";
    AddHelpLine('d', "d", "delta parameter", std::to_string(delta_));
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'd':
        if (std::is_floating_point<WeightT_>::value)
          delta_ = static_cast<WeightT_>(atof(opt_arg));
        else
          delta_ = static_cast<WeightT_>(atol(opt_arg));
        break;
      default: CLApp::HandleArg(opt, opt_arg);
    }
  }

  WeightT_ delta() const { return delta_; }
};



class CLConvert : public CLBase {
  std::string out_filename_ = "";
  bool out_weighted_ = false;
  bool out_el_ = false;
  bool out_sg_ = false;

 public:
  CLConvert(int argc, char** argv, std::string name)
      : CLBase(argc, argv, name) {
    get_args_ += "e:b:w";
    AddHelpLine('b', "file", "output serialized graph to file");
    AddHelpLine('e', "file", "output edge list to file");
    AddHelpLine('w', "file", "make output weighted");
  }

  void HandleArg(signed char opt, char* opt_arg) override {
    switch (opt) {
      case 'b': out_sg_ = true; out_filename_ = std::string(opt_arg);   break;
      case 'e': out_el_ = true; out_filename_ = std::string(opt_arg);   break;
      case 'w': out_weighted_ = true;                                   break;
      default: CLBase::HandleArg(opt, opt_arg);
    }
  }

  std::string out_filename() const { return out_filename_; }
  bool out_weighted() const { return out_weighted_; }
  bool out_el() const { return out_el_; }
  bool out_sg() const { return out_sg_; }
};

#endif  // COMMAND_LINE_H_
