#pragma once
#include "common.h"
#include "Sequence.h"


class SeedFinder {
    virtual void CreateIndexFromGenome(const Sequence& genome) = 0;
    virtual void GetSeedsWithRead(const Sequence& read, vector<match>& seeds) = 0;
};