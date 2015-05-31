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

HighErrorRateRegions herRegions;
int numConnectable = 0;
int numConnectableBefore = 0;

int outputCounter = 0;
string GetFilename() {
    char buffer[100];
    sprintf(buffer, "seeds_alignments/%07d.png", outputCounter++);
    return string(buffer);
}

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

        
        if (Utility::IsSignificant(al)) {
            alList.push_back(al);
        }
        
    }
    
    for (auto it = alList.begin(); it != alList.end(); ++it) {
        bool doRealign = false;
        //Visualization vis("");
        for (auto it2 = alList.begin(); it2 != alList.end(); ++it2) {
            if (it == it2) continue;
            if (Utility::AreConnectable(*it, *it2)) {
                doRealign = true;
                numConnectableBefore++;
                //vis.AddAlignment(*it);
                //vis.CreateVisualization("", 2, 500, 500);
                //vis.AddAlignment(*it2, cv::Scalar(255, 0, 0));
                //vis.CreateVisualization("", 2, 500, 500);
                //vis.CreateVisualization(GetFilename(), 2, 500, 500);
                //herRegions.AddRegion(*it, *it2);
                break;
            }
        }
        if (doRealign) {
            Alignment newAlignment;
            dwHighTolerance.ComputeAlignment(genome, read, it->GetStartPos(), newAlignment);
            newAlignment.ComputeTrace();
            //vis.AddAlignment(newAlignment, cv::Scalar(0, 255, 0));
            //vis.CreateVisualization(GetFilename(), 2, 500, 500);
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
        longestAlignment = max(longestAlignment, al.GetLengthOnB());
        op->AddAlignment(al, read, dir);
    }
    
    if (longestAlignment < read.GetData().length() * 0.8) {
        for (auto it = alList.begin(); it != alList.end(); ++it) {
            for (auto it2 = alList.begin(); it2 != alList.end(); ++it2) {
                if (it == it2) continue;
                if (Utility::AreConnectable(*it, *it2)) {
                    herRegions.AddRegion(*it, *it2);
                    numConnectable++;
                    break;
                }
            }
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
    dw.SetAligningParameters(0.7, 50, {0.25, 0.25, 0.25, 0.25});
    dwHighTolerance.SetAligningParameters(0.6, 50, {0.25, 0.25, 0.25, 0.25});

    // find seeds
    SeedFinder_hashmap_2bit_1hashmap sf(15);
    sf.CreateIndexFromGenome(genome);
        
    double seedsSum = 0;
    double expandSum = 0;
    
    while (fastq >> read) {
        vector<Match> forwardMatches;
        vector<Match> backwardMatches;
        
        sf.GetSeedsWithRead(read, forwardMatches, backwardMatches);

        vector<Alignment> forwardAlignments, backwardAlignments;
        int longestForward = ExpandSeeds(genome, read, forwardMatches, forwardAlignments, EDIR_FORWARD);

        int longestBackward = ExpandSeeds(reverseGenome, read, backwardMatches, backwardAlignments, EDIR_BACKWARD);
        
        Visualization vis1;
        for (auto a : forwardAlignments) {
            vis1.AddAlignment(a);
        }
        vis1.AddSeeds(forwardMatches);
        vis1.CreateVisualization(GetFilename(), 2, 1500, 500, genome.GetData().length(), read.GetData().length());
        
        Visualization vis2;
        for (auto a : backwardAlignments) {
            vis2.AddAlignment(a);
        }
        vis2.AddSeeds(backwardMatches);
        vis2.CreateVisualization(GetFilename(), 2, 1500, 500, genome.GetData().length(), read.GetData().length());
    }
    cout << "num connectable: " << numConnectable << endl;
    cout << "num connectable before: " << numConnectableBefore << endl;
    herRegions.OutputHERRegions(3);
    delete op;
}

