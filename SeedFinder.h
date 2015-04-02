#pragma once
#include "common.h"
#include "Sequence.h"

struct match {
    int genomePos, readPos, length;
};

class SeedFinder {
    virtual void CreateIndexFromGenome(const Sequence& genome) = 0;
    virtual void GetSeedsWithRead(const Sequence& read, vector<match>& seeds) = 0;
};