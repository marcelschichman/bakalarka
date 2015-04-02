#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdio>
#define FOR(II, NN) for (int II = 0; II < NN; II++)

struct match {
    int genomePos, readPos, length;
};