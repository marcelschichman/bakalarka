#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Utility.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
#include "BlasrOutput.h"
#include "SeedFinder_hashmap_2bit.h"
#include "SeedFinder_hashmap_2bit_1hashmap.h"
#include <list>
#include <fstream>
#include <set>
using namespace std;

Output *op;
DalignWrapper dw;
DalignWrapper dwHighTolerance;
Visualization *viz;



int ExpandSeeds(Sequence &genome, Sequence &read, vector<Match> &matches, vector<Alignment> &alignments, Direction dir) {
    AlignedPairsSet alignedPairs;
    
    int longestAlignment = 0;
    int alignmentsCounter = 0;
//    cout << "matches " << matches.size() << endl;
//    VisualizationData vd;
    list<Alignment> alList;
    for (Match &match : matches) {
        if (alignedPairs.IsAligned(match)) {
            continue;
        }
        Alignment al;
        dw.ComputeAlignment(genome, read, match, al);

        vector<pair<int, int>> pairs;
        al.ComputeTrace();
        al.GetAlignedPairs(pairs);
        alignedPairs.MarkAsAligned(pairs);

        longestAlignment = max(longestAlignment, al.GetLengthOnB());

//        if (al.GetLengthOnB() >= read.GetData().length() - 20) {
//            op->AddAlignment(al, read, dir);
//        }
//        if (al.GetLengthOnB() > read.GetData().length() * 0.1) {
//            viz->AddAlignment(al);
//        }
        
        if (Utility::IsSignificant(al)) {
            alList.push_back(al);
        }
        
    }
    
    for (auto it = alList.begin(); it != alList.end(); ++it) {
        bool doRealign = false;
//        Visualization vis(read.GetId());
        for (auto it2 = alList.begin(); it2 != alList.end(); ++it2) {
            if (it == it2) continue;
            if (Utility::AreConnectable(*it, *it2)) {
                doRealign = true;
//                vis.AddAlignment(*it);
//                vis.CreateVisualization("", 2, 500, 500);
//                vis.AddAlignment(*it2, cv::Scalar(255, 0, 0));
//                vis.CreateVisualization("", 2, 500, 500);
                break;
            }
        }
        if (doRealign) {
            Alignment newAlignment;
            dwHighTolerance.ComputeAlignment(genome, read, it->GetStartPos(), newAlignment);
            newAlignment.ComputeTrace();
//            vis.AddAlignment(newAlignment, cv::Scalar(0, 255, 0));
//            vis.CreateVisualization("", 2, 500, 500);
            *it = newAlignment;
        }
        
        vector<pair<int, int>> pairs;
        it->GetAlignedPairs(pairs);
        AlignedPairsSet alignedPairsSingleAl;
        alignedPairsSingleAl.MarkAsAligned(pairs);
        
        list<Alignment>::iterator erasable;
        for (auto it2 = alList.begin(); it2 != alList.end(); ) {
            erasable = it2++;
            if (it == erasable) continue;
            if (alignedPairsSingleAl.IsAligned(erasable->GetStartPos()) || alignedPairsSingleAl.IsAligned(erasable->GetEndPos())) {
                alList.erase(erasable);
            }
        }
    }
    alignments.assign(alList.begin(), alList.end());
    
    for (auto &al : alignments) {
        op->AddAlignment(al, read, dir);
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
    dwHighTolerance.SetAligningParameters(0.6, 50, {0.25, 0.25, 0.25, 0.25});

    // find seeds
    SeedFinder_hashmap_2bit_1hashmap sf(20);
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
//        viz = new Visualization();
//        
//        if (forwardMatches.size() > backwardMatches.size()) {
//            int longestForward = ExpandSeeds(genome, read, forwardMatches, EDIR_FORWARD);
//        } else {
//            int longestBackward = ExpandSeeds(reverseGenome, read, backwardMatches, EDIR_BACKWARD);
//        }
//        
//        viz->AddSeeds(forwardMatches.size() > backwardMatches.size() ? forwardMatches : backwardMatches);
//        
//        viz->CreateVisualization("", 1, 1800, 300, genome.GetData().length(), read.GetData().length());
//        delete viz;
        vector<Alignment> forwardAlignments, backwardAlignments;
        int longestForward = ExpandSeeds(genome, read, forwardMatches, forwardAlignments, EDIR_FORWARD);

        int longestBackward = ExpandSeeds(reverseGenome, read, backwardMatches, backwardAlignments, EDIR_BACKWARD);
        
//        for (Alignment &al : forwardAlignments.size() > backwardAlignments.size() ? forwardAlignments : backwardAlignments) {
//            Visualization vis(read.GetId());
//            vis.AddAlignment(al);
//            vis.AddSeeds(forwardAlignments.size() > backwardAlignments.size() ? forwardMatches : backwardMatches);
//            vis.CreateVisualization("", 2, 1800, 300, genome.GetData().length(), read.GetData().length());
//        }
        
        //expandSum += Utility::GetTimerResult();
        //cout << read.GetId() << " " << max(longestForward, longestBackward) << endl;
        //Utility::verbalResult(to_string(read.GetData().length()) + " long read ");
    }
//    clog << "seeds: " << seedsSum << endl;
//    clog << "expand: " << expandSum << endl;
    delete op;
}

