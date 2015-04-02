#pragma once
#include "common.h"
#include "Sequence.h"

struct match {
    int genomePos, readPos, length;
    operator pair<int, int>() {
        return make_pair(genomePos, readPos);
    }
};


class SeedFinder {
    virtual void CreateIndexFromGenome(const Sequence& genome) = 0;
    virtual void GetSeedsWithRead(const Sequence& read, vector<match>& seeds) = 0;
};