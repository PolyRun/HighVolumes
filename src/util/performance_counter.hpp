// Performance_Counter Utility
//
// Build a recursive structure corresponding to your code.
// Make sure you count for the correct function choice.
//
// PC_Stack: keeps track of the function stack, and how often the current funciton was called.
//           maps function choice (void*) to PC_Cost_Wrapper, which holds the cost function
// 
// Register functions:
// pc_stack().add((void*)someFunction, new PC_Cost_Wrapper<someFunction_cost_f>(someFunction_cost,"someFunction"));
//
// Call cost function:
//   pc_stack().reset();
//   {
//      PC_Frame<someCost_cost_f> frame((void*)someFunction); // provide function or variable holding function choice
//      frame.costf()(cost_function_inputs);
//   }
//   pc_stack().print();


#ifndef HEADER_PERFORMANCE_COUNTER_HPP
#define HEADER_PERFORMANCE_COUNTER_HPP

#include <cassert>
#include <map>
#include <iostream>
#include <deque>
#include "table.hpp"

class PC_Cost_Wrapper_Base {
   void virtual placeholder() {};// just to make things polymorphic ;)
};
template <class CF_t>
class PC_Cost_Wrapper : public PC_Cost_Wrapper_Base {
public:
   PC_Cost_Wrapper(CF_t costf, const std::string& name) : costf(costf), name(name) {}
   CF_t costf;
   std::string name;
};

class PC_Stack {
public:
   PC_Stack() : table(Table({"Stack","split","calls","flops","bytes","Flops","Bytes","tag"})) {}
   void add(void* f, PC_Cost_Wrapper_Base* wrapper) {
      functions[f] = wrapper;
   }
   PC_Cost_Wrapper_Base* get(void* f) {
      auto it = functions.find(f);
      if(it == functions.end()) {
         return NULL;
      } else {
         return it->second;
      }
   };
   void log(const size_t _flops, const size_t _bytes, const std::string& tag) {
      std::cout << "  " << getPath() 
	        << " (" << _flops << " " << _bytes << ") "
	        << " (" << _flops*current_reps << " " << _bytes*current_reps << ") "
		<< tag <<  "\n";
      
      table.addRow({
		    {"Stack",getPath()},
		    {"bytes",std::to_string(_bytes)},
		    {"Bytes",std::to_string(_bytes*current_reps)},
		    {"flops",std::to_string(_flops)},
		    {"Flops",std::to_string(_flops*current_reps)},
		    {"tag",tag}});

      flops_ += _flops*current_reps;
      bytes_ += _bytes*current_reps;
   }
   void reset() {
      flops_ = 0;
      bytes_ = 0;
      current_reps = 1;
      reps.clear();
      frames.clear();
   }
   void print() {
      std::cout << "PC_Stack: flops: " << flops_ << ", bytes: " << bytes_ << "\n";
      table.print();
   }
   size_t flops() {return flops_;}
   size_t bytes() {return bytes_;}
   void pushFrame(const std::string& name, size_t reps_) {
      reps.push_back(current_reps);
      current_reps *= reps_;
      frames.push_back(name);
      std::cout << "> " << getPath() << " " << reps_ << " " << current_reps <<  "\n";
      table.addRow({
		    {"Stack",getPath()},
		    {"split",std::to_string(reps_)},
		    {"calls",std::to_string(current_reps)}});
   }
   void popFrame() {
      std::cout << "< " << getPath() <<  "\n";
      current_reps = reps.back();
      reps.pop_back();
      frames.pop_back();
   }
   std::string getPath() {
      std::string ret = "";
      for(auto s : frames) {
         ret+= "/" + s;
      }
      return ret;
   }
private:
   std::map<void*, PC_Cost_Wrapper_Base*> functions;
   size_t flops_ = 0;
   size_t bytes_ = 0;
   std::deque<size_t> reps;
   std::deque<std::string> frames;
   size_t current_reps = 1;
   Table table;
};

PC_Stack& pc_stack();

class PC_Frame_Base {
public:
   PC_Frame_Base() {}
   PC_Frame_Base(const std::string& name, size_t reps) {
      pc_stack().pushFrame(name, reps);
   }
   ~PC_Frame_Base() {
      pc_stack().popFrame();
   }
};
template <class CF_t>
class PC_Frame : public PC_Frame_Base {
public:
   PC_Frame(void* f, size_t reps=1) {
      PC_Cost_Wrapper_Base* b = pc_stack().get(f);
      assert(b && "must find cost_wrapper for f! - did you register all functions?");
      wrapper = dynamic_cast<PC_Cost_Wrapper<CF_t>*>(b);
      assert(wrapper && "cost_wrapper must have correct cost function type!");
      pc_stack().pushFrame(name(), reps);
   }
   CF_t costf() {return wrapper->costf;}
   const std::string& name() {return wrapper->name;}
private:
   PC_Cost_Wrapper<CF_t>* wrapper = NULL;
};

#endif // HEADER_PERFORMANCE_COUNTER_HPP
