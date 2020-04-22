// Header only for Timer

#include <chrono>

extern "C" {
#include "tsc_x86.h"
}

#ifndef HEADER_TIMER_HPP
#define HEADER_TIMER_HPP


class Timer_generic{
public:   
    virtual void start() = 0;
    virtual void stop() = 0;
    virtual double get_time() = 0;
};


class Timer : public Timer_generic {
public:
    Timer() {};
    void start() {t1_=t0_=std::chrono::high_resolution_clock::now();}
    void stop() {t1_=std::chrono::high_resolution_clock::now();}
    double get_time(){ return millisecs(); }
    double millisecs() {return std::chrono::duration_cast<std::chrono::duration<double,std::milli>>(t1_-t0_).count();}
    double microsecs() {return std::chrono::duration_cast<std::chrono::duration<double,std::micro>>(t1_-t0_).count();}
    double seconds() {return std::chrono::duration_cast<std::chrono::duration<double>>(t1_-t0_).count();}
private:
    std::chrono::high_resolution_clock::time_point t0_, t1_;
};


class Tsc : public Timer_generic {
public:
    Tsc() {};
    void start() {
        t0 = start_tsc();
        tdiff = 0;
    }
    void stop(){
        tdiff = stop_tsc(t0);
    }
    double get_time(){
        return (double) tdiff;
    }

private:
    myInt64 t0, tdiff;
};

#endif // HEADER_TIMER_HPP
