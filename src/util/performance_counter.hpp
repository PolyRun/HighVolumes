
#ifndef HEADER_PERFORMANCE_COUNTER_HPP
#define HEADER_PERFORMANCE_COUNTER_HPP

#include <cassert>
#include <map>
#include <iostream>

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
   void log(const int _flops, const int _bytes) {
      std::cout << "Log " << _flops << " " << _bytes << "\n";
      flops_ += _flops;
      bytes_ += _bytes;
   }
   int flops() {return flops_;}
   int bytes() {return bytes_;}
private:
   std::map<void*, PC_Cost_Wrapper_Base*> functions;
   int flops_ = 0;
   int bytes_ = 0;
};

PC_Stack& pc_stack();

class PC_Frame_Base {};
template <class CF_t>
class PC_Frame : public PC_Frame_Base {
public:
   PC_Frame(void* f) {
      PC_Cost_Wrapper_Base* b = pc_stack().get(f);
      wrapper = dynamic_cast<PC_Cost_Wrapper<CF_t>*>(b);
      assert(wrapper && "must find cost_wrapper for f!");
      std::cout << "Frame start: " << name() << "\n";
   }
   ~PC_Frame() {
      std::cout << "Frame end: " << name() << "\n";
   }
   CF_t costf() {return wrapper->costf;}
   const std::string& name() {return wrapper->name;}
private:
   PC_Cost_Wrapper<CF_t>* wrapper = NULL;
};

#endif // HEADER_PERFORMANCE_COUNTER_HPP
