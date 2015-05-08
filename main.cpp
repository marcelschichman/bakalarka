#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Utility.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
#include "BlasrOutput.h"
#include "SeedFinder_hashmap_2bit.h"
#include <forward_list>
#include <fstream>
#include <set>
using namespace std;

Output *op;
DalignWrapper dw;
Visualization *viz;



int ExpandSeeds(Sequence &genome, Sequence &read, vector<Match> &matches, Direction dir) {
    set<pair<int, int>> alignedPairs;
    
    int longestAlignment = 0;
    int alignmentsCounter = 0;
//    cout << "matches " << matches.size() << endl;
//    VisualizationData vd;
    for (Match &match : matches) {
        if (alignedPairs.find(make_pair(match.genomePos+8, match.readPos+8)) != alignedPairs.end()) {
            continue;
        }
        Alignment al;
        dw.ComputeAlignment(genome, read, match, al);

        vector<pair<int, int>> pairs;
        al.ComputeTrace();
        al.GetAlignedPairs(pairs);
        alignedPairs.insert(pairs.begin(), pairs.end());

        longestAlignment = max(longestAlignment, al.GetLengthOnB());
//        pair<int, int> posA = al.GetPosOnA();
//        pair<int, int> posB = al.GetPosOnB();
//        cout << al.GetLengthOnB() << " posA: " << posA.first << " " << posA.second << " posb: " << posB.first << " " << posB.second << endl;

//        if (al.GetLengthOnB() >= read.GetData().length() - 20) {
//            op->AddAlignment(al, read, dir);
//            break;
//        }
//        if (al.GetLengthOnB() > read.GetData().length() * 0.1) {
            viz->AddAlignment(al);
//        }
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
    dw.SetAligningParameters(0.7, 50, {0.25, 0.25, 0.25, 0.25});

    // find seeds
    SeedFinder_hashmap_2bit sf(20);
    sf.CreateIndexFromGenome(genome);
        
    double seedsSum = 0;
    double expandSum = 0;
    
    while (fastq >> read) {
        
        //Utility::startTiming();
        vector<Match> forwardMatches;
        vector<Match> backwardMatches;
        //Utility::StartTiming();
        sf.GetSeedsWithRead(read, forwardMatches, backwardMatches);
        //seedsSum += Utility::GetTimerResult();
//        
//        auto &matches = forwardMatches.size() > backwardMatches.size() ? forwardMatches : backwardMatches;
//        Utility::Filter(forwardMatches, 100, 0.3);
//        Utility::Filter(backwardMatches, 100, 0.3);
//        for_each(matches.begin(), matches.end(), [](Match m) { cout << m.score << " "; });
//        cout << endl;
        
        //Utility::StartTiming();
        //Utility::Filter();
        viz = new Visualization();
        
        if (forwardMatches.size() > backwardMatches.size()) {
            int longestForward = ExpandSeeds(genome, read, forwardMatches, EDIR_FORWARD);
        } else {
            int longestBackward = ExpandSeeds(reverseGenome, read, backwardMatches, EDIR_BACKWARD);
        }
        
        viz->AddSeeds(forwardMatches.size() > backwardMatches.size() ? forwardMatches : backwardMatches);
        
        viz->CreateVisualization("", 1, 1800, 300, genome.GetData().length(), read.GetData().length());
        delete viz;
        //expandSum += Utility::GetTimerResult();
        //cout << read.GetId() << " " << max(longestForward, longestBackward) << endl;
        //Utility::verbalResult(to_string(read.GetData().length()) + " long read ");
    }
//    clog << "seeds: " << seedsSum << endl;
//    clog << "expand: " << expandSum << endl;
    delete op;
}

