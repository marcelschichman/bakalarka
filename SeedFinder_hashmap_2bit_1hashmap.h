#pragma once
#include "common.h"
#include "SeedFinder.h"
#include <unordered_map>
using namespace std;

class SeedFinder_hashmap_2bit_1hashmap : public SeedFinder {
public:
    SeedFinder_hashmap_2bit_1hashmap(int _length) :
    length(_length) {
    }

    virtual void CreateIndexFromGenome(Sequence& genome);
    virtual void GetSeedsWithRead(Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds);
protected:
    int length;
    unordered_map<long long, vector<pair<Direction, int>>> kMerMap;
    
    void ExpandPairs(vector<pair<int, int>>& kMerPairs, vector<Match>& seeds);
    static bool comparePairsByDiagonal(const pair<int, int> &left, const pair<int, int> &right);
};
