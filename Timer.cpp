#include "Timer.h"

vector<chrono::time_point<chrono::system_clock>> Timer::tpts;

double Timer::getTimerResult() {
    double result = ((chrono::duration<double>)(chrono::system_clock::now() - tpts.back())).count();
    tpts.pop_back();
    return result;
}

void Timer::startTiming() {
    tpts.push_back(chrono::system_clock::now());
}

void Timer::verbalResult(const string& name) {
    cout << name << " took " << getTimerResult() << "s" << endl;
}
