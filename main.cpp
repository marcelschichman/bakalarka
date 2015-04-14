#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Timer.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
#include <forward_list>
#include <set>
using namespace std;

int THE_APP(int argc, char** argv) {
    string genomeFilename = "genome.fasta";
    string readsFilename = "pacbio_10kb.fastq";
    string outputFilename = "output.txt";
    
    if (argc == 4) {
        genomeFilename = argv[1];
        readsFilename = argv[2];
        outputFilename = argv[3];
    }
    
    // get sequences
    Sequence a, b;
    FASTA fasta(genomeFilename);
    FASTQ fastq(readsFilename);
    fasta >> a;
    
    SAMOutput soutput(outputFilename);

    // prepare wrapper
    DalignWrapper dw;
    dw.SetAligningParameters(0.70, 5, {0.25, 0.25, 0.25, 0.25});

    // find seeds
    SeedFinder_hashmap sf(30);
    sf.CreateIndexFromGenome(a);
        
    int readCounter = 0;
    while (fastq >> b) {
        if (b.data.length() < 1000) {
            continue;
        }
        cout << readCounter++ << " ";
        
        vector<Match> matches;
        sf.GetSeedsWithRead(b, matches);
        
        set<pair<int, int>> alignedPairs;
        
        int alignmentsCounter = 0;
        for (Match &match : matches) {
            if (alignedPairs.find(make_pair(match.genomePos+15, match.readPos+15)) != alignedPairs.end()) {
                continue;
            }
            Alignment al;
            dw.ComputeAlignment(a, b, match, al);
            cout << al.GetLengthOnB() << " " << b.data.length() << " " << (matches.front().genomePos - matches.front().readPos) << endl;
            
            vector<pair<int, int>> pairs;
            al.ComputeTrace();
            al.GetAlignedPairs(pairs);
            alignedPairs.insert(pairs.begin(), pairs.end());
            
            if (al.GetLengthOnB() == b.data.length()) {
                soutput.AddAlignment(b, al);
            }
        }
    }
}

