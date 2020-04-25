
#ifndef HEADER_PERFORMANCE_COUNTER_HPP
#define HEADER_PERFORMANCE_COUNTER_HPP

#include <cassert>
#include <map>
#include <iostream>
#include <deque>

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
   PC_Stack() {}
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
      flops_ += _flops*current_reps;
      bytes_ += _bytes*current_reps;
   }
   void reset() {
      flops_ = 0;
      bytes_ = 0;
   }
   void print() {
      std::cout << "PC_Stack: flops: " << flops_ << ", bytes: " << bytes_ << "\n";
   }
   size_t flops() {return flops_;}
   size_t bytes() {return bytes_;}
   void pushFrame(const std::string& name, size_t reps_) {
      reps.push_back(current_reps);
      current_reps *= reps_;
      frames.push_back(name);
      std::cout << "> " << getPath() << " " << reps_ << " " << current_reps <<  "\n";
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
      wrapper = dynamic_cast<PC_Cost_Wrapper<CF_t>*>(b);
      assert(wrapper && "must find cost_wrapper for f!");
      pc_stack().pushFrame(name(), reps);
   }
   CF_t costf() {return wrapper->costf;}
   const std::string& name() {return wrapper->name;}
private:
   PC_Cost_Wrapper<CF_t>* wrapper = NULL;
};

#endif // HEADER_PERFORMANCE_COUNTER_HPP
