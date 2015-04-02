#pragma once
#include "common.h"
#include <chrono>
using namespace std;
namespace Timer {
    vector<chrono::time_point<chrono::system_clock>> tpts;

    void startTiming() {
        tpts.push_back(chrono::system_clock::now());
    }

    double getTimerResult() {
        double result = ((chrono::duration<double>)(chrono::system_clock::now() - tpts.back())).count();
        tpts.pop_back();
        return result;
    }
};
