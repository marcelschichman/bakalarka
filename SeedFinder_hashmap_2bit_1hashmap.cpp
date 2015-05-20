#include "SeedFinder_hashmap_2bit_1hashmap.h"

void SeedFinder_hashmap_2bit_1hashmap::CreateIndexFromGenome(Sequence& genome) {
    long long seed = 0;
    long long reversedSeed = 0;
    long long mask = ((long long)1 << (2 * length)) - 1;
    char *data = genome.ToDalignFromat();
    int genomeLength = (int)genome.GetData().length();
    kMerMap.reserve(genomeLength);
    FOR(i, genomeLength) {
        seed = ((seed << 2) + data[i]) & mask;
        reversedSeed = ((reversedSeed << 2) + 3 - data[genomeLength - i - 1]) & mask;
        
        if (i >= length - 1) {
            kMerMap[seed].push_back(make_pair(Direction::EDIR_FORWARD, i - length + 1));
            kMerMap[reversedSeed].push_back(make_pair(Direction::EDIR_BACKWARD, i - length + 1));
        }
    }
}

bool SeedFinder_hashmap_2bit_1hashmap::comparePairsByDiagonal(const pair<int, int> &left, const pair<int, int> &right) {
    int leftDiagonal = left.second - left.first;
    int rightDiagonal = right.second - right.first;
    return (leftDiagonal == rightDiagonal) ? (left.first < right.first) : (leftDiagonal < rightDiagonal);
}

void SeedFinder_hashmap_2bit_1hashmap::GetSeedsWithRead(Sequence& read, vector<Match>& forwardSeeds, vector<Match>& reverseSeeds) {
    vector<pair<int, int>> kMerPairs;
    vector<pair<int, int>> kMerPairsReversed;
    
    long long seed = 0;
    long long mask = ((long long)1 << (2 * length)) - 1;
    char *data = read.ToDalignFromat();

    FOR(i, (int)read.GetData().length()) {
        seed = ((seed << 2) + data[i]) & mask;
        
        if (i >= length) {
            auto positions = kMerMap.find(seed);
            if (positions != kMerMap.end()) {
                for (auto &position : positions->second) {
                    if (position.first == Direction::EDIR_FORWARD) {
                        kMerPairs.push_back(make_pair(position.second, i - length + 1));
                    } else {
                        kMerPairsReversed.push_back(make_pair(position.second, i - length + 1));
                    }
                }
            }
        }
    }
    
    if (kMerPairs.size() == 0 && kMerPairsReversed.size() == 0) {
        length = (length / 2) * 2;
        long long firstHalf = 0;
        long long secondHalf = 0;
        long long mask = ((long long)1 << length) - 1;
        FOR(i, (int)read.GetData().length() - (length / 2) - 1) {
            firstHalf = ((firstHalf << 2) + data[i]) & mask;
            secondHalf = ((secondHalf << 2) + data[i + (length / 2) + 1]) & mask;
            
            if (i >= length / 2) {
                seed = secondHalf + (firstHalf << length);
                auto positions = kMerMap.find(seed);
                if (positions != kMerMap.end()) {
                    //cout << "zufale seedy: " << read.GetId() << endl;
                    for (auto &position : positions->second) {
                        //cout << position.first << " " << position.second << " " << i - (length / 2) + 1  << endl;
                        if (position.first == Direction::EDIR_FORWARD) {
                            kMerPairs.push_back(make_pair(position.second, i - (length / 2) + 1));
                        } else {
                            kMerPairsReversed.push_back(make_pair(position.second, i - (length / 2) + 1));
                        }
                    }
                }
            }
        }
    }

    ExpandPairs(kMerPairs, forwardSeeds);
    ExpandPairs(kMerPairsReversed, reverseSeeds);
}

void SeedFinder_hashmap_2bit_1hashmap::ExpandPairs(vector<pair<int, int> >& kMerPairs, vector<Match>& seeds) {
    sort(kMerPairs.begin(), kMerPairs.end(), SeedFinder_hashmap_2bit_1hashmap::comparePairsByDiagonal);

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
