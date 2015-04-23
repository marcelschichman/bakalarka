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

    virtual void CreateIndexFromGenome(const Sequence& genome);
    virtual void GetSeedsWithRead(const Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds);
protected:
    int length;
    unordered_map<string, vector<int>> kMerMap;
    unordered_map<string, vector<int>> kMerMapReversed;
    
    void ExpandPairs(vector<pair<int, int>>& kMerPairs, vector<Match>& seeds);
};
