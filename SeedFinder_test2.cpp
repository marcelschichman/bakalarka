#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Utility.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
#include "SeedFinder_hashmap_2bit.h"
#include "SeedFinder_hashmap_2bit_1hashmap.h"
#include <forward_list>
#include <set>
#include <cstdlib>
#include <unordered_map>
using namespace std;

int SEED_FINDER_TEST2(int argc, char** argv) {
    
    FASTA fasta("genome.fasta");
    FASTQ fastq("1000_long_pbsim_reads.fastq");
    
    Sequence genome, read;
    fasta >> genome;
    unordered_map<long long, vector<pair<Direction, int>>> kMerMap;
    int length = 20;
    {
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
    
    vector<pair<int, int>> numSeeds;
    while (fastq >> read) {
        pair<int, int> newValue;
        
        ///// normalne seedy
        long long seed = 0;
        long long mask = ((long long)1 << (2 * length)) - 1;
        char *data = read.ToDalignFromat();

        FOR(i, (int)read.GetData().length()) {
            seed = ((seed << 2) + data[i]) & mask;

            if (i >= length - 1) {
                auto positions = kMerMap.find(seed);
                if (positions != kMerMap.end()) {
                    for (auto &position : positions->second) {
                        newValue.first++;
                    }
                }
            }
        }
        
        ///// zufale seedy
        long long firstHalf = 0;
        long long secondHalf = 0;
        mask = ((long long)1 << length) - 1;
        FOR(i, (int)read.GetData().length() - (length / 2) - 1) {
            firstHalf = ((firstHalf << 2) + data[i]) & mask;
            secondHalf = ((secondHalf << 2) + data[i + (length / 2) + 1]) & mask;
            
            if (i >= length / 2 - 1) {
                seed = secondHalf + (firstHalf << length);
                auto positions = kMerMap.find(seed);
                if (positions != kMerMap.end()) {
                    for (auto &position : positions->second) {
                        newValue.second++;
                    }
                }
            }
        }
        numSeeds.push_back(newValue);
    }
    
    ofstream os("seedfinder_comparison.txt");
    int both = 0;
    int none = 0;
    int only_normal = 0;
    int only_desperate = 0;
    for (auto x : numSeeds) {
        os << x.first << " " << x.second << endl;
        
        if (x.first > 0) {
            if (x.second > 0)
                both++;
            else
                only_normal++;
        }
        else {
            if (x.second > 0)
                only_desperate++;
            else
                none++;
        }
    }
    cout << "genome length: " << genome.GetData().length() << endl;
    cout << "total: " << (both + none + only_desperate + only_normal) << endl;
    cout << "both: " << both << endl;
    cout << "none: " << none << endl;
    cout << "only normal: " << only_normal << endl;
    cout << "only desperate: " << only_desperate << endl;
}
