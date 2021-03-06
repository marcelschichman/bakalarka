#pragma once
#include "common.h"
#include "SeedFinder.h"
#include <unordered_map>
using namespace std;

class SeedFinder_hashmap : public SeedFinder {
public:
    SeedFinder_hashmap(int _length) :
    length(_length) {
    }

    virtual void CreateIndexFromGenome(Sequence& genome);
    virtual void GetSeedsWithRead(Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds);
protected:
    int length;
    unordered_map<string, vector<int>> kMerMap;
    unordered_map<string, vector<int>> kMerMapReversed;
    
    void ExpandPairs(vector<pair<int, int>>& kMerPairs, vector<Match>& seeds);
    static bool comparePairsByDiagonal(pair<int, int> left, pair<int, int> right);
};
