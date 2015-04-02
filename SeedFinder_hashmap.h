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
    virtual void GetSeedsWithRead(const Sequence& read, vector<match>& seeds);
private:
    int length;
    unordered_map<string, vector<int>> kMerMap;
};
