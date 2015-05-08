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
    virtual void CreateIndexFromGenome(Sequence& genome) = 0;
    virtual void GetSeedsWithRead(Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds) = 0;
};