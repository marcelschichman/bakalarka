#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Timer.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
#include "BlasrOutput.h"
#include "MatchFilter.h"
#include <forward_list>
#include <set>
using namespace std;

Output *op;
DalignWrapper dw;



int ExpandSeed(Sequence &genome, Sequence &read, vector<Match> &matches, Direction dir) {
    set<pair<int, int>> alignedPairs;
    
    int longestAlignment = 0;
    int alignmentsCounter = 0;
    cout << "matches " << matches.size() << endl;
    for (Match &match : matches) {
        if (alignedPairs.find(make_pair(match.genomePos+4, match.readPos+4)) != alignedPairs.end()) {
            continue;
        }
        Alignment al;
        dw.ComputeAlignment(genome, read, match, al);

        vector<pair<int, int>> pairs;
        al.ComputeTrace();
        al.GetAlignedPairs(pairs);
        alignedPairs.insert(pairs.begin(), pairs.end());
        int numEqual = 0;

        longestAlignment = max(longestAlignment, al.GetLengthOnB());
        pair<int, int> pos = al.GetPosOnB();
        cout << al.GetLengthOnB() << " " << pos.first << " " << pos.second << endl;

        if (al.GetLengthOnB() >= read.GetData().length() - 20) {
            op->AddAlignment(al, read, dir);
            break;
        }
    }
    return longestAlignment;
}

int THE_APP(int argc, char** argv) {
    string genomeFilename = "genome.fasta";
    string readsFilename = "fajny_read.fastq";
    string outputFilename = "output.txt";
    
    if (argc == 4) {
        genomeFilename = argv[1];
        readsFilename = argv[2];
        outputFilename = argv[3];
    }
    
    // get sequences
    Sequence genome, read;
    FASTA fasta(genomeFilename);
    FASTQ fastq(readsFilename);
    fasta >> genome;
    Sequence reverseGenome (genome, true);
    
    op = new BlasrOutput(outputFilename, genome);
    op->Init();

    // prepare wrapper
    dw.SetAligningParameters(0.6, 5, {0.25, 0.25, 0.25, 0.25});

    // find seeds
    SeedFinder_hashmap sf(15);
    sf.CreateIndexFromGenome(genome);
        
    double seedsSum = 0;
    double expandSum = 0;
    
    while (fastq >> read) {
        //Timer::startTiming();
        vector<Match> forwardMatches;
        vector<Match> backwardMatches;
        Timer::startTiming();
        sf.GetSeedsWithRead(read, forwardMatches, backwardMatches);
        seedsSum += Timer::getTimerResult();
//        
//        auto &matches = forwardMatches.size() > backwardMatches.size() ? forwardMatches : backwardMatches;
        MatchFilter mfForward(forwardMatches);
        mfForward.Sort(100);
        mfForward.Filter(0.5);
        MatchFilter mfBackward(backwardMatches);
        mfBackward.Sort(100);
        mfBackward.Filter(0.5);
//        for_each(matches.begin(), matches.end(), [](Match m) { cout << m.score << " "; });
//        cout << endl;
        
        Timer::startTiming();
        int longestForward = ExpandSeed(genome, read, forwardMatches, EDIR_FORWARD);
        cout << endl;
        //int longestBackward = ExpandSeed(reverseGenome, read, backwardMatches, EDIR_BACKWARD);
        expandSum += Timer::getTimerResult();
        //cout << read.GetId() << " " << max(longestForward, longestBackward) << endl;
        //Timer::verbalResult(to_string(read.GetData().length()) + " long read ");
    }
    clog << "seeds: " << seedsSum << endl;
    clog << "expand: " << expandSum << endl;
    delete op;
}

