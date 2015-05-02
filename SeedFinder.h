#pragma once
#include "common.h"
#include "Sequence.h"

struct Match {
    int genomePos, readPos, length, score;
    operator pair<int, int>() {
        return make_pair(genomePos, readPos);
    }
};

class SeedFinder {
    virtual void CreateIndexFromGenome(const Sequence& genome) = 0;
    virtual void GetSeedsWithRead(const Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds) = 0;
};