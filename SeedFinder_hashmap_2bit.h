#pragma once
#include "common.h"
#include "SeedFinder.h"
#include <unordered_map>
using namespace std;

class SeedFinder_hashmap_2bit : public SeedFinder {
public:
    SeedFinder_hashmap_2bit(int _length) :
    length(_length) {
    }

    virtual void CreateIndexFromGenome(Sequence& genome);
    virtual void GetSeedsWithRead(Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds);
protected:
    int length;
    unordered_map<long long, vector<int>> kMerMap;
    unordered_map<long long, vector<int>> kMerMapReversed;
    
    void ExpandPairs(vector<pair<int, int>>& kMerPairs, vector<Match>& seeds);
    static bool comparePairsByDiagonal(const pair<int, int>& left, const pair<int, int>& right);
};
