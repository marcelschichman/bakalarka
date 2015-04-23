#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Timer.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
#include "BlasrOutput.h"
#include <forward_list>
#include <set>
using namespace std;

Output *op;
DalignWrapper dw;



int ExpandSeed(Sequence &genome, Sequence &read, vector<Match> &matches, Direction dir) {
    set<pair<int, int>> alignedPairs;
    
    int longestAlignment = 0;
    int alignmentsCounter = 0;
    for (Match &match : matches) {
        if (alignedPairs.find(make_pair(match.genomePos+15, match.readPos+15)) != alignedPairs.end()) {
            continue;
        }
        Alignment al;
        dw.ComputeAlignment(genome, read, match, al);

        vector<pair<int, int>> pairs;
        al.ComputeTrace();
        al.GetAlignedPairs(pairs);
        alignedPairs.insert(pairs.begin(), pairs.end());

        longestAlignment = max(longestAlignment, al.GetLengthOnB());

        if (al.GetLengthOnB() >= read.GetData().length() - 20) {
            op->AddAlignment(al, read, dir);
            break;
        }
    }
    return longestAlignment;
}

int THE_APP(int argc, char** argv) {
    string genomeFilename = "genome.fasta";
    string readsFilename = "100_3000.fastq";
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
    dw.SetAligningParameters(0.70, 5, {0.25, 0.25, 0.25, 0.25});

    // find seeds
    SeedFinder_hashmap sf(15);
    sf.CreateIndexFromGenome(genome);
        
    while (fastq >> read) {
        Timer::startTiming();
        vector<Match> forwardMatches;
        vector<Match> backwardMatches;
        sf.GetSeedsWithRead(read, forwardMatches, backwardMatches);
        
        int longestForward = ExpandSeed(genome, read, forwardMatches, EDIR_FORWARD);
        int longestBackward = ExpandSeed(reverseGenome, read, backwardMatches, EDIR_BACKWARD);
        //cout << max(longestForward, longestBackward) << endl;
        Timer::verbalResult(to_string(read.GetData().length()) + " long read ");
    }
    delete op;
}

