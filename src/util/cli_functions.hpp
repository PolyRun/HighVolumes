// Header only for CLI functions
//
// Your situation
//  you have multiple functions that perform the same job, you need to choose one
//  you would like to do the choice via CLI input
//
// How it works:
//  hooks up to CLI, makes sure it filters the right inputs
//  then, takes input from CLI, sets correct function pointers to the functions of choice
//
// Implementation:
//  CLIFunctions -> configure and commit (preParse) before parse of cli
//  after parse, read arguments (postParse), set function pointers
//
// Recommendation:
//  Extend from CLIFunctions, make config in constructor

#include <set>

#include "cli.hpp"


#ifndef HEADER_CLI_FUNCTIONS_HPP
#define HEADER_CLI_FUNCTIONS_HPP

class CLIF_OptionBase {
public:
   CLIF_OptionBase(const signed char opt, const std::string &name, const std::string &val)
	   : name_(name), opt_(opt), val_(val) {}
   // called by CLIFunctions in postParse
   // sets function pointer to option chosen by cli
   void virtual postParse(CLI &cli) {}

   const signed char opt_; // opt char for cli input
   const std::string name_; // name of parameter under opt
   const std::string val_; // default value for parameter
};

template <class F_t>
class CLIF_Option : public CLIF_OptionBase {
public:
   CLIF_Option(F_t* var, const signed char opt, const std::string &name, const std::string &val,
	       const std::map<std::string, F_t> &m)
   : var_(var), CLIF_OptionBase(opt,name,val), fmap(m) {}
   
   void virtual postParse(CLI &cli) {
      std::string key = cli.parameters(opt_).get(name_,"");
      const auto it = fmap.find(key);
      if(it==fmap.end()) {
         std::cout << "Error: bad function name for: " << opt_ << " " << name_ << "\n";
	 std::cout << "       viable options:\n";
	 for(const auto it : fmap) {
	    std::cout << "           " << it.first << "\n";
	 }
	 std::exit(0);
      }
      *var_ = fmap.at(cli.parameters(opt_).get(name_,""));
   }
   const std::map<std::string, F_t> fmap; // map holding all the options
private:
   F_t* var_; // global variable that will hold final choice
};

class CLIFunctions {
public:
   // initialized this before cli is parsed
   CLIFunctions(CLI &cli) : cli_(cli) {
      // initialize function lists
      // set cli parameters accordingly
   }

   // call this after cli was parsed
   void postParse() {
      // read cli -> set function ptrs
   
      for(auto o : options_) {
         o->postParse(cli_);
      }
   }

   // register function ptr
   // char + name
   // list of functions + their names
   void add(CLIF_OptionBase *o) {
      auto it = params_.find(o->opt_);
      if(it==params_.end()) {
         // check if opt is still available:
         if(cli_.isUsed(o->opt_)) {
	    std::cout << "Error: CLIFunctions.add: opt " << o->opt_ << " not free!\n";
	    std::exit(0);
	 }
	 params_[o->opt_] = CLIParameters();
      }
      
      // register the parameter
      CLIParameters &p = params_[o->opt_];
      if(p.isUsed(o->name_)) {
         std::cout << "Error: CLIFunction.add opt " << o->opt_ << " parameter " << o->name_ << " not free!\n";
	 std::exit(0);
      }
      p.set(o->name_,o->val_);
      
      options_.insert(o);
   }
   
   // call before cli.parse, to commit all registered functions
   void preParse() {
      for(auto it : params_) {
	 // register the CLIParameters for opt
	 cli_.addParameters(it.first,it.second, "Function configuration");
      }
   }

   CLIF_OptionBase* getOption(const std::string &name) {
      for(auto it : options_) {
         if(it->name_ == name) {return it;}
      }
      return NULL;
   }

private:
   CLI &cli_;
   std::map<signed char,CLIParameters> params_; // assembled parameters for opt
   std::set<CLIF_OptionBase*> options_;
};

#endif // HEADER_CLI_FUNCTIONS_HPP

