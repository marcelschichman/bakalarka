#include "SeedFinder_hashmap_2bit.h"

void SeedFinder_hashmap_2bit::CreateIndexFromGenome(Sequence& genome) {
    long long seed = 0;
    long long reversedSeed = 0;
    long long mask = ((long long)1 << (2 * length)) - 1;
    char *data = genome.ToDalignFromat();
    int genomeLength = (int)genome.GetData().length();
    FOR(i, genomeLength) {
        seed = ((seed << 2) + data[i]) & mask;
        reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;
        
        if (i >= length - 1) {
            kMerMap[seed].push_back(i - length + 1);
            kMerMapReversed[reversedSeed].push_back(i - length + 1);
        }
    }
}

bool SeedFinder_hashmap_2bit::comparePairsByDiagonal(const pair<int, int> &left, const pair<int, int> &right) {
    int leftDiagonal = left.second - left.first;
    int rightDiagonal = right.second - right.first;
    return (leftDiagonal == rightDiagonal) ? (left.first < right.first) : (leftDiagonal < rightDiagonal);
}

void SeedFinder_hashmap_2bit::GetSeedsWithRead(Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds) {
    vector<pair<int, int>> kMerPairs;
    vector<pair<int, int>> kMerPairsReversed;
    
    long long seed = 0;
    long long mask = ((long long)1 << (2 * length)) - 1;
    char *data = read.ToDalignFromat();

    FOR(i, (int)read.GetData().length()) {
        seed = ((seed << 2) + data[i]) & mask;
        
        if (i >= length) {
            auto positions = kMerMap.find(seed);
            auto positionsReversed = kMerMapReversed.find(seed);
            if (positions != kMerMap.end()) {
                for (int position : positions->second) {
                    kMerPairs.push_back(make_pair(position, i - length + 1));
                }
            }
            if (positionsReversed != kMerMapReversed.end()) {
                for (int position : positionsReversed->second) {
                    kMerPairsReversed.push_back(make_pair(position, i - length + 1));
                }
            }
        }
    }

    ExpandPairs(kMerPairs, forwardSeeds);
    ExpandPairs(kMerPairsReversed, reverseSeeds);
}

void SeedFinder_hashmap_2bit::ExpandPairs(vector<pair<int, int> >& kMerPairs, vector<Match>& seeds) {
    sort(kMerPairs.begin(), kMerPairs.end(), SeedFinder_hashmap_2bit::comparePairsByDiagonal);

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
