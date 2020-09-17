// Header only for CLI

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <getopt.h>
#include <algorithm>
#include <iomanip>

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

// -------------------- How to use CLI_LONG:
// CLI_LONG cli(argc,argv,"peterem");
// cli.addFlag('x',"xopt","Some test flag");
// cli.addFlag('y',"yopt","Some test flag");
// cli.addOption('z',"zopt","defaultv","Some test option");
// cli.addOption('Z',"z-opt-two","defaultv","Some test option");
// cli.parse();
// 
// std::cout << "flag x: " << cli.flag('x') << "\n";
// std::cout << "flag y: " << cli.flag('y') << "\n";
// std::cout << "option z: " << cli.option('z') << "\n";
// std::cout << "option Z: " << cli.option('Z') << "\n";
//
// ---- disabled parameters (dict) because ugly.


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

class CLI_LONG {
public:
   CLI_LONG(int argc, char** argv, std::string name = "")
   : argc_(argc), argv_(argv), name_(name) {
      addFlag('h',"help","print this help message");
   }
   
   // call after all flags, options, params declared
   bool parse() {
      signed char c; // character
      extern char *optarg; // getopt - optional string
      int option_index = 0;

      long_options_.push_back({0,0,0,0});
      
      while ((c = getopt_long(argc_, argv_, parse_.c_str(),
		              long_options_.data(), &option_index)) != -1) {
	 // std::cout << "getopt " << c << "\n";
	 if(c=='?' || !handleArg(c,optarg)) {
	    usage();
	    std::exit(0);
	 }
      }
      
     
      /* Print any remaining command line arguments (not options). */
      if (optind < argc_) {
	 std::cout << "Unexpected arguments:\n";
         while (optind < argc_)
            std::cout << "  " << argv_[optind++] << "\n";
         usage();
	 std::exit(0);
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
         auto it = flags_.find(opt2long_[opt]);
         if(it!=flags_.end()) {
            it->second = !(it->second);
            return true;
         }
      }
      {// try find option:
         auto it = option_.find(opt2long_[opt]);
         if(it!=option_.end()) {
            it->second = std::string(opt_arg);
            return true;
         }
      }
      
      std::cout << "arg not found: " << opt << "\n";
      return false; // nothing found.
   }
   
   // print usage to screen
   void usage() {
      std::cout << "Usage: " << name_ << std::endl;

      for(std::map<std::string, signed char>::iterator it = long2opt_.begin(); it!=long2opt_.end(); it++) {
	 std::cout << " -" << it->second << " --";
	 std::cout << std::left << std::setw(15) << it->first;
	 std::cout << " " << desc_[it->first];
	 
	 auto itf = flags_.find(it->first);
	 if(itf!=flags_.end()) {
	    std::cout << " (flag)";
	 }
	 
	 std::cout << std::endl;
      }
   }
   
   void checkOpt(signed char opt) {
      auto it = opt2long_.find(opt);
      if(it!=opt2long_.end()) {
         std::cout << "Error: -" << opt << " already exists!\n";
         std::exit(0);
      }
   }
   void checkOptLong(const std::string &optlong) {
      auto it = long2opt_.find(optlong);
      if(it!=long2opt_.end()) {
         std::cout << "Error: --" << optlong << " already exists!\n";
         std::exit(0);
      }
   }

   void addFlag(signed char opt, const std::string &optlong, std::string desc) {
      checkOpt(opt);
      checkOptLong(optlong);
      desc_.insert(std::pair<std::string,std::string>(optlong,desc));
      flags_.insert(std::pair<std::string,bool>(optlong,false));
      parse_ += opt;
      
      opt2long_.insert(std::pair<signed char,std::string>(opt,optlong));
      long2opt_.insert(std::pair<std::string,signed char>(optlong,opt));
      
      // did not just take optlong, because it may disappear. opt2long_[opt] should persist.
      long_options_.push_back({opt2long_[opt].c_str(),no_argument,0,opt});
   }

   bool flag(const std::string &optlong) {
      auto it = flags_.find(optlong);
      if(it!=flags_.end()) {
         return it->second;
      } else {
         std::cout << "Error: flag " << optlong << " never declared!\n";
         return false;
      }
   }

   void addOption(signed char opt, const std::string &optlong, std::string def, std::string desc) {
      checkOpt(opt);
      checkOptLong(optlong);
         
      desc += " (default: " + def + ")";
      desc_.insert(std::pair<std::string,std::string>(optlong,desc));
      option_.insert(std::pair<std::string,std::string>(optlong,def));
      parse_ += opt;
      parse_ += ":";
   
      opt2long_.insert(std::pair<signed char,std::string>(opt,optlong));
      long2opt_.insert(std::pair<std::string,signed char>(optlong,opt));
      
      // did not just take optlong, because it may disappear. opt2long_[opt] should persist.
      long_options_.push_back({opt2long_[opt].c_str(),required_argument,0,opt});
   }

   std::string option(const std::string &optlong) {
      auto it = option_.find(optlong);
      if(it!=option_.end()) {
         return it->second;
      } else {
         std::cout << "Error: option " << optlong << " never declared!\n";
         return "";
      }
   }

   bool isUsed(signed char opt) {return opt2long_.find(opt)!=opt2long_.end();}
   bool isUsedLong(const std::string &optlong) {return long2opt_.find(optlong)!=long2opt_.end();}

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
   // arguments for main
   const int argc_;
   char** argv_;
   
   // App name
   const std::string name_;
   
   ///  // string listing the short opt characters
   std::string parse_ = "";
   // array of structs for long options
   std::vector<struct option> long_options_;
   
   std::map<std::string, std::string> desc_; // desc for usage display
   std::map<signed char, std::string> opt2long_; // long option name
   std::map<std::string, signed char> long2opt_; // long option name
   std::map<std::string, bool> flags_; // bool for flags, only exists if flag exists
   std::map<std::string, std::string> option_; /// string for a option, exists if option exists
};



#endif // HEADER_CLI_HPP

