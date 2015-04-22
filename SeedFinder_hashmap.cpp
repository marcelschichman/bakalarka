#include "SeedFinder_hashmap.h"

void SeedFinder_hashmap::CreateIndexFromGenome(const Sequence& genome) {
    FOR(i, (int)genome.GetData().length() - length + 1) {
        kMerMap[genome.GetData().substr(i, length)].push_back(i);
    }
}

bool comparePairsByDiagonal(pair<int, int> left, pair<int, int> right) {
    int leftDiagonal = left.second - left.first;
    int rightDiagonal = right.second - right.first;
    return (leftDiagonal == rightDiagonal) ? (left.first < right.first) : (leftDiagonal < rightDiagonal);
}

void SeedFinder_hashmap::GetSeedsWithRead(const Sequence& read, vector<Match>& seeds) {
    vector<pair<int, int>> kMerPairs;

    FOR(i, (int)read.GetData().length() - length + 1) {
        string kMer = read.GetData().substr(i, length);
        auto positions = kMerMap.find(kMer);
        if (positions != kMerMap.end()) {
            for (int position : positions->second) {
                kMerPairs.push_back(make_pair(position, i));
            }
        }
    }

    sort(kMerPairs.begin(), kMerPairs.end(), comparePairsByDiagonal);

    pair<int, int> previous = make_pair(-5, -5);
    for (auto ppair : kMerPairs) {
        if (previous.first == ppair.first - 1 && previous.second == ppair.second - 1) {
            seeds.back().length++;
        } else {
            Match newMatch = {ppair.first, ppair.second, length};
            seeds.push_back(newMatch);
        }
        previous = ppair;
    }
}

