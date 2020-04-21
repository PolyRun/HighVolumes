// Header only for CLI

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <getopt.h>
#include <algorithm>

// -------------------- Idea:
//
// Takes input from argc, argv
// can be configured to interpret it as flags, options (strings) or parameters (dictionary)
// can set default values
//
// -------------------- How to use:
//
// CLI cli(argc,argv,"AppName");
//
// cli.addFlag('v',false,"verbose"); // cli.flag('v')
// cli.addOption('n',"default","size"); // cli.option('n')
// 
// CLIParameters clip;
// clip.set("p1","param_value");
// cli.addParameters('p',clip,"somee extra parameters");
// std::cout << cli.parameters('p') << std::endl; // print all
// cli.parameters('p').get("p1","default");


#ifndef HEADER_CLI_HPP
#define HEADER_CLI_HPP


// Basically implements a dictionary
class CLIParameters {
public:
   CLIParameters() {}
   CLIParameters(std::string in_) {parse(in_);}
   void parse(std::string in_) {
      std::vector<std::string> parts = split(in_, ',');
      for(std::string m : parts) {
         std::vector<std::string> pp = split(m,'=');
	 //assert(pp.size()==2 && "parameter list must have correct form!");
	 if(! (pp.size()==2 && "parameter list must have correct form!")) {
	    std::cout << "Parameter list had bad form (expect comma seperated list of var=val):\n" << in_ << "\n";
	    std::exit(0);
	 }
	 params_[pp[0]] = pp[1];
      }
   }
   std::string get(const std::string &name, const std::string &value) const {
      auto it = params_.find(name);
      if(it==params_.end()) {
	 return value;
      } else {
         return it->second;
      }
   }
   void set(const std::string &name, const std::string &value) {
      params_[name] = value;
   }

   std::string tostring() const {
      std::string ret;
      for(auto it = params_.begin(); it!=params_.end(); it++) {
         ret+= "    " + it->first + "=" + it->second + "\n";
      }
      
      return ret;
   }

   bool isUsed(const std::string &name) {return params_.find(name)!=params_.end();}
private:
   std::map<std::string,std::string> params_;
   
   std::vector<std::string> split(std::string &s, char separator =',') {
       std::vector<std::string> res;
       std::istringstream f(s);
       std::string tmp;
       while (getline(f, tmp, separator)) {
           res.push_back(tmp);
       }
       return res;
   }
};

class CLI {
public:
   CLI(int argc, char** argv, std::string name = "")
   : argc_(argc), argv_(argv), name_(name) {
      addFlag('h',false,"print this help message");
   }
   
   // call after all flags, options, params declared
   bool parse() {
      signed char c_opt; // character
      extern char *optarg; // getopt - optional string
      while ((c_opt = getopt(argc_, argv_, parse_.c_str())) != -1) {
         if(!handleArg(c_opt,optarg)) {return false;}
      }
      return true;
   }
   
   // internal: handles input from cl, assigns to flags, options, params
   bool virtual handleArg(signed char opt, char* opt_arg) {
      if(opt=='h') {
	 usage();
	 std::exit(0);
      }
      {// try find flag:
         auto it = flags_.find(opt);
	 if(it!=flags_.end()) {
	    it->second = !(it->second);
	    return true;
	 }
      }
      {// try find option:
         auto it = option_.find(opt);
	 if(it!=option_.end()) {
	    it->second = std::string(opt_arg);
	    return true;
	 }
      }
      {// try find params:
         auto it = params_.find(opt);
	 if(it!=params_.end()) {
	    it->second.parse(std::string(opt_arg));
	    return true;
	 }
      }
     
      return false; // nothing found.
   }
   
   // print usage to screen
   void usage() {
      std::cout << "Usage: " << name_ << std::endl;
      for(std::map<signed char, std::string>::iterator it = desc_.begin(); it!=desc_.end(); it++) {
         std::cout << " " << it->first << " " << it->second << std::endl;
      }
   }

   void addFlag(signed char opt, bool def, std::string desc) {
      auto it = desc_.find(opt);
      if(it!=desc_.end()) {
         std::cout << "Error: addFlag " << opt << " already exists!\n";
	 std::exit(0);
      } else {
         desc_.insert(std::pair<signed char,std::string>(opt,desc));
         flags_.insert(std::pair<signed char,bool>(opt,def));
	 parse_ += opt;
      }
   }

   bool flag(signed char opt) {
      auto it = flags_.find(opt);
      if(it!=flags_.end()) {
	 return it->second;
      } else {
         std::cout << "Error: flag " << opt << " never declared!\n";
	 return false;
      }
   }

   void addOption(signed char opt, std::string def, std::string desc) {
      auto it = desc_.find(opt);
      if(it!=desc_.end()) {
         std::cout << "Error: addOption " << opt << " already exists!\n";
	 std::exit(0);
      } else {
	 desc += "\n    default: " + def;
         desc_.insert(std::pair<signed char,std::string>(opt,desc));
         option_.insert(std::pair<signed char,std::string>(opt,def));
	 parse_ += opt;
	 parse_ += ":";
      }
   }

   std::string option(signed char opt) {
      auto it = option_.find(opt);
      if(it!=option_.end()) {
	 return it->second;
      } else {
         std::cout << "Error: option " << opt << " never declared!\n";
	 return "";
      }
   }

   void addParameters(signed char opt, const CLIParameters &def, std::string desc) {
      auto it = desc_.find(opt);
      if(it!=desc_.end()) {
         std::cout << "Error: addParameters " << opt << " already exists!\n";
	 std::exit(0);
      } else {
	 desc += "\n" + def.tostring();
         desc_.insert(std::pair<signed char,std::string>(opt,desc));
         params_.insert(std::pair<signed char,CLIParameters>(opt,def));
	 parse_ += opt;
	 parse_ += ":";
      }
   }

   CLIParameters& parameters(signed char opt) {
      auto it = params_.find(opt);
      if(it!=params_.end()) {
	 return it->second;
      } else {
         std::cout << "Error: parameters " << opt << " never declared!\n";
	 return nullParams;
      }
   }

   bool isUsed(signed char opt) {return desc_.find(opt)!=desc_.end();}

    std::string getPath(){
        std::string s(argv_[0]);
        return s;
    }

    std::string getPathFromExec(){
        std::string path = getPath();
        std::string res = "";
        reverse(path.begin(), path.end());
        size_t pos = path.find('/');
        // the executable is not in the current directory
        if(pos != std::string::npos){
            reverse(path.begin(), path.end());
            res = path.substr(0, path.length() - pos);
        }
        return res;
    }
    
protected:
   const int argc_;
   char** argv_;
   const std::string name_;

   std::string parse_ = "";
   
   std::map<signed char, std::string> desc_;
   std::map<signed char, bool> flags_;
   std::map<signed char, std::string> option_;
   std::map<signed char, CLIParameters> params_;
private:
   CLIParameters nullParams; // if none set
};

#endif // HEADER_CLI_HPP

