#include "SeedFinder_hashmap.h"

void SeedFinder_hashmap::CreateIndexFromGenome(const Sequence& genome) {
    Sequence reversedGenome (genome, true);
    FOR(i, (int)genome.GetData().length() - length + 1) {
        kMerMap[genome.GetData().substr(i, length)].push_back(i);
        kMerMapReversed[reversedGenome.GetData().substr(i, length)].push_back(i);
    }
}

bool comparePairsByDiagonal(pair<int, int> left, pair<int, int> right) {
    int leftDiagonal = left.second - left.first;
    int rightDiagonal = right.second - right.first;
    return (leftDiagonal == rightDiagonal) ? (left.first < right.first) : (leftDiagonal < rightDiagonal);
}

void SeedFinder_hashmap::GetSeedsWithRead(const Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds) {
    vector<pair<int, int>> kMerPairs;
    vector<pair<int, int>> kMerPairsReversed;

    FOR(i, (int)read.GetData().length() - length + 1) {
        string kMer = read.GetData().substr(i, length);
        auto positions = kMerMap.find(kMer);
        auto positionsReversed = kMerMapReversed.find(kMer);
        if (positions != kMerMap.end()) {
            for (int position : positions->second) {
                kMerPairs.push_back(make_pair(position, i));
            }
        }
        if (positionsReversed != kMerMapReversed.end()) {
            for (int position : positionsReversed->second) {
                kMerPairsReversed.push_back(make_pair(position, i));
            }
        }
    }

    ExpandPairs(kMerPairs, forwardSeeds);
    ExpandPairs(kMerPairsReversed, reverseSeeds);
}

void SeedFinder_hashmap::ExpandPairs(vector<pair<int, int> >& kMerPairs, vector<Match>& seeds) {
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
