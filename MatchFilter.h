#pragma once
#include "common.h"
#include "SeedFinder.h"

class MatchFilter {
public:
    MatchFilter(vector<Match> &_matches);
    void Sort(int nearDist);
    void Filter(double remove);
private:
    vector<Match> &matches;
    int Diagonal(Match &match);
};
