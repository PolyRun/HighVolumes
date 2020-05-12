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
#include <sstream>

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
	       const std::map<std::string, std::pair<F_t,std::string>> &m)
   : var_(var), CLIF_OptionBase(opt,name,val), fmap(m) {}
   
   void virtual postParse(CLI &cli) {
      std::string key = cli.parameters(opt_).get(name_,"");
      const auto it = fmap.find(key);
      if(it==fmap.end()) {
          std::cout << "Error: bad choice " << key << " for: " << opt_ << " " << name_ << "\n";
	 std::cout << "       viable options:\n";
	 for(const auto it : fmap) {
	    std::cout << "           " << it.first << " - " << it.second.second << "\n";
	 }
	 std::exit(0);
      }
      *var_ = fmap.at(cli.parameters(opt_).get(name_,"")).first;
   }
   std::map<std::string, std::pair<F_t,std::string>> fmap; // map holding all the options
   // maps name to {value, description}
private:
   F_t* var_; // global variable that will hold final choice
};


template <class F_t, class G_t>
class CLIF_DoubleOption : public CLIF_OptionBase {
public:
   CLIF_DoubleOption(F_t* varF,G_t* varG, const signed char opt, const std::string &name, const std::string &val,
	       const std::map<std::string, std::pair<std::pair<F_t,G_t>,std::string>> &m)
   : varF_(varF), varG_(varG), CLIF_OptionBase(opt,name,val), fmap(m) {}
   
   void virtual postParse(CLI &cli) {
      std::string key = cli.parameters(opt_).get(name_,"");
      const auto it = fmap.find(key);
      if(it==fmap.end()) {
          std::cout << "Error: bad choice " << key << " for: " << opt_ << " " << name_ << "\n";
	 std::cout << "       viable options:\n";
	 for(const auto it : fmap) {
	    std::cout << "           " << it.first << " - " << it.second.second << "\n";
	 }
	 std::exit(0);
      }
      *varF_ = fmap.at(cli.parameters(opt_).get(name_,"")).first.first;
      *varG_ = fmap.at(cli.parameters(opt_).get(name_,"")).first.second;
   }
   std::map<std::string, std::pair<std::pair<F_t,G_t>,std::string>> fmap; // map holding all the options
   // maps name to {value, description}
private:
   F_t* varF_; // global variable that will hold final choice
   G_t* varG_; // global variable that will hold final choice
};


template <class F_t, class G_t, class H_t>
class CLIF_TrippleOption : public CLIF_OptionBase {
public:
   CLIF_TrippleOption(F_t* varF,G_t* varG,H_t* varH, const signed char opt, const std::string &name, const std::string &val,
	       const std::map<std::string, std::pair<std::pair<F_t,std::pair<G_t,H_t>>,std::string>> &m)
   : varF_(varF), varG_(varG), varH_(varH), CLIF_OptionBase(opt,name,val), fmap(m) {}
   
   void virtual postParse(CLI &cli) {
      std::string key = cli.parameters(opt_).get(name_,"");
      const auto it = fmap.find(key);
      if(it==fmap.end()) {
          std::cout << "Error: bad choice " << key << " for: " << opt_ << " " << name_ << "\n";
	 std::cout << "       viable options:\n";
	 for(const auto it : fmap) {
	    std::cout << "           " << it.first << " - " << it.second.second << "\n";
	 }
	 std::exit(0);
      }
      *varF_ = fmap.at(cli.parameters(opt_).get(name_,"")).first.first;
      *varG_ = fmap.at(cli.parameters(opt_).get(name_,"")).first.second.first;
      *varH_ = fmap.at(cli.parameters(opt_).get(name_,"")).first.second.second;
   }
   std::map<std::string, std::pair<std::pair<F_t,std::pair<G_t,H_t>>,std::string>> fmap; // map holding all the options
   // maps name to {value, description}
private:
   F_t* varF_; // global variable that will hold final choice
   G_t* varG_; // global variable that will hold final choice
   H_t* varH_; // global variable that will hold final choice
};




template <class T>
class CLIF_OptionNumber : public CLIF_OptionBase {
public:
   CLIF_OptionNumber(T* var, const signed char opt, const std::string &name, const std::string &val,
	       T minVal, T maxVal)
   : var_(var), CLIF_OptionBase(opt,name,val), minVal_(minVal), maxVal_(maxVal) {}
   
   void virtual postParse(CLI &cli) {
      std::string in = cli.parameters(opt_).get(name_,"");
      
      T value;
      std::stringstream convert(in);
      convert >> value;
      
      if(value > maxVal_ or value < minVal_) {
         std::cout << "Error: parameter out of limits for: " << opt_ << " " << name_ << "\n";
	 std::cout << "       viable range: " << minVal_ << " to " << maxVal_ << "\n";
	 std::exit(0);
      }
      *var_ = value;
   }
   T minVal_, maxVal_;
private:
   T* var_; // global variable that will hold final choice
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
   
   // claim opt code for a specific description
   void claimOpt(const signed char opt, const std::string &desc) {
      auto it = params_.find(opt);
      if(it==params_.end()) {
         // check if opt is still available:
         if(cli_.isUsed(opt)) {
	    std::cout << "Error: CLIFunctions.claimOpt: opt " << opt << " not free!\n";
	    std::exit(0);
	 }
	 params_[opt] = CLIParameters();
	 desc_[opt] = desc;
      } else {
         std::cout << "Error: CLIFunctions.claimOpt: opt " << opt << " not free!\n";
	 std::exit(0);
      }
   }

   // register function ptr
   // char + name
   // list of functions + their names
   void add(CLIF_OptionBase *o) {
      auto it = params_.find(o->opt_);
      if(it==params_.end()) {
         // check if opt is still available:
         claimOpt(o->opt_, "Function configuration");
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
	 cli_.addParameters(it.first,it.second, desc_[it.first]);
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
   std::map<signed char,std::string> desc_; 
   std::set<CLIF_OptionBase*> options_;
};

#endif // HEADER_CLI_FUNCTIONS_HPP

