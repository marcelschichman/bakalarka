#pragma once
#include "common.h"
#include "Sequence.h"

struct Match {
    int genomePos, readPos, length;
    operator pair<int, int>() {
        return make_pair(genomePos, readPos);
    }
};


class SeedFinder {
    virtual void CreateIndexFromGenome(const Sequence& genome) = 0;
    virtual void GetSeedsWithRead(const Sequence& read, vector<Match>& seeds) = 0;
};