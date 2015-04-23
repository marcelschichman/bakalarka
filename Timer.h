#pragma once
#include "common.h"
#include <chrono>

//#ifndef MY_TIMER_INCLUDE
//#define MY_TIMER_INCLUDE
using namespace std;
class Timer {
    static vector<chrono::time_point<chrono::system_clock>> tpts;

public:
    static void startTiming();
    static double getTimerResult();
    static double verbalResult(const string& name);
};