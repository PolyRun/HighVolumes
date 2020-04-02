// Header only for Timer

#include <chrono>

#ifndef HEADER_TIMER_HPP
#define HEADER_TIMER_HPP

class Timer {
public:
   Timer(){};
   void start() {t1_=t0_=std::chrono::high_resolution_clock::now();}
   void stop() {t1_=std::chrono::high_resolution_clock::now();}
   double millisecs() {return std::chrono::duration_cast<std::chrono::duration<double,std::milli>>(t1_-t0_).count();}
   double microsecs() {return std::chrono::duration_cast<std::chrono::duration<double,std::micro>>(t1_-t0_).count();}
   double seconds() {return std::chrono::duration_cast<std::chrono::duration<double>>(t1_-t0_).count();}
private:
   std::chrono::high_resolution_clock::time_point t0_, t1_;
};

#endif // HEADER_TIMER_HPP
