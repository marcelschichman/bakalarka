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

double Timer::verbalResult(const string& name) {
    double result = getTimerResult();
    cout << name << " took " << result << "s" << endl;
    return result;
}
