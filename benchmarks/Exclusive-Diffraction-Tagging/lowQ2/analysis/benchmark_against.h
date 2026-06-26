#pragma once
#include "/opt/benchmarks/common_bench/include/common_bench/benchmark.h"

namespace common_bench {
  //want to create benchmark tests for each variable
  //either from scratch or using a previous test as target value *some constant
  
  //{"name", "title", "description", "quantity",
  //"target", "value", "result"}

  bool does_file_exist(const std::string& name){
    ifstream f(name.c_str());
    return f.good();
  }
  
  class test_against {

  public:
    
  test_against(const std::string& oldfile,double tolerance=1.01 ):_tolerance(tolerance){

      std::ifstream f(oldfile);
      _json=nlohmann::json::parse(f);
      cout<<_json["tests"][0]["name"]<<endl;
 
    }

    
    void add_gt_test(const std::string& testName, double testValue){
      find_test(testName);
      //my target should be the old test value
      auto target = get_value()/_tolerance;
      _json["tests"][_curr_test]["target"]=fmt::format("{}", target);
     _json["tests"][_curr_test]["result"]= testValue > target ? "pass" : "fail";
      std::cout<<"Result of test "<<testName<<_json["tests"][_curr_test]["result"]<<" "<<testValue<<" > "<<target<<std::endl;
      _json["tests"][_curr_test]["value"]=testValue;

    }
    
    void add_lt_test(const std::string& testName, double testValue){
      find_test(testName);
      //my target should be the old test value
      auto target = get_value()*_tolerance;
      _json["tests"][_curr_test]["target"]=fmt::format("{}", target);
      _json["tests"][_curr_test]["result"]= testValue < target ? "pass" : "fail";
      std::cout<<"Result of test "<<testName<<_json["tests"][_curr_test]["result"]<<" "<<testValue<<" < "<<target<<std::endl;
     _json["tests"][_curr_test]["value"]=testValue;
    }
    
    std::string get_name(){return _json["tests"][_curr_test]["name"];}
    std::string get_title(){return _json["tests"][_curr_test]["title"];}
    std::string get_description(){return _json["tests"][_curr_test]["description"];}
    std::string get_quantity(){return _json["tests"][_curr_test]["quantity"];}
    std::string get_target_string(){return _json["tests"][_curr_test]["target"];}
    double get_target(){return std::strtod(get_target_string().data(), nullptr);}
    double get_value(){return _json["tests"][_curr_test]["value"];}
    double get_weight(){return _json["tests"][_curr_test]["weight"];}
    std::string get_result(){return _json["tests"][_curr_test]["result"];}

    void find_test(const std::string& tname){
      _curr_test=0;
      for(auto& test : _json["tests"]){
	if(test["name"]==tname) return;
	++_curr_test;
      }
      std::cerr<<"findTest not found exiting ...."<<tname<<std::endl;exit(0);
    }
    void write_tests(const std::string& newFile){
      std::cout << fmt::format("Writing test data to {}\n", newFile);
      std::ofstream output_file(newFile);
      output_file << std::setw(4) << _json << "\n";

    }
  private :
    
    nlohmann::json _json;
    uint _curr_test=0; //index of current test, found via find_test
    double _tolerance = 1.01; //increase target by this factor
  };

  void bench_template(const std::string& fname, uint ntests=1){
    common_bench::Test test1{{
	{"name", "name"},
	  {"title", "title"},
	    {"description", "description"},
	      {"quantity", "quantity"},
		{"target", "target"}}};
    std::vector<common_bench::Test> tests(ntests,test1);
    common_bench::write_test(tests, fname);
  }

  
}
