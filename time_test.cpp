#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Utility.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
#include "SeedFinder_hashmap_2bit.h"
#include <forward_list>
#include <set>
#include <cstdlib>
using namespace std;

vector<int> groups = {10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1400, 1600, 1800, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 20000, 1 << 30};
vector<vector<double>> timesAlignment (groups.size());
vector<vector<double>> timesTrace (groups.size());
vector<vector<double>> timesPointPairs (groups.size());
vector<vector<double>> timesFindSeeds (groups.size());
vector<vector<double>> timesFindInAligned (groups.size());

int GetGroup(int x);
int ExpandReads(Sequence &genome, Sequence &read, vector<Match> &matches, DalignWrapper &dw);

#define NUM_READS 100000

int TIME_TEST(int argc, char** argv) {
    //header
    cout << "*** Execution time test ***" << endl;
    cout << "This tests tests how much the length of the alignment influences the execution time of function calls." << endl;
    cout << "One time operations: " << endl;
    
    string genomeFilename = "genome.fasta";
    string readsFilename = "pb_1000_plus.fastq";
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
    
    Utility::StartTiming();
    fasta >> a;
    Utility::VerbalResult("Reading genome");
    cout << "Genome length: " << a.GetData().length() << endl;
    Sequence ra (a, true);
    
    //SAMOutput soutput(outputFilename);

    // prepare wrapper
    Utility::StartTiming();
    DalignWrapper dw;
    dw.SetAligningParameters(0.70, 50, {0.25, 0.25, 0.25, 0.25});
    Utility::VerbalResult("Preparing DALIGN");

    // find seeds
    Utility::StartTiming();
    SeedFinder_hashmap_2bit sf(20);
    sf.CreateIndexFromGenome(a);
    Utility::VerbalResult("Preparing seed finder");
        
    int readCounter = 0;
    int progressBarPos = 0;
    
    while (fastq >> b) {
        // get seed finder time
        Utility::StartTiming();
        vector<Match> forwardMatches;
        vector<Match> backwardMatches;
        sf.GetSeedsWithRead(b, forwardMatches, backwardMatches);
        timesFindSeeds[GetGroup(b.GetData().size())].push_back(Utility::GetTimerResult());   
        
//        Utility::Filter(forwardMatches, backwardMatches, 300, 0.5);
        
        ExpandReads(a, b, forwardMatches, dw);
        ExpandReads(ra, b, backwardMatches, dw);
        
        if (readCounter % 1000 == 0) 
            clog << readCounter << endl;
        if (++readCounter >= NUM_READS)
            break;
    }
    
    cout << "Per read operations: " << endl;
    cout << "group\tnum\tseeds\taligned" << endl;
    for (int i = 0; i < groups.size(); i++) {
        if (timesFindSeeds[i].size() == 0) {
            continue;
        }
        double findSeedsAverage = 0;
        for (double x : timesFindSeeds[i]) {
            findSeedsAverage += x;
        }
        cout << groups[i] << "\t";
        cout << timesFindSeeds[i].size() << "\t";
        cout << (findSeedsAverage / timesFindSeeds[i].size()) << "\t";
        
        if (timesFindInAligned[i].size()) {
            double findInPairedAverage = 0;
            for (double x : timesFindInAligned[i]) {
                findInPairedAverage += x;
            }
            cout << (findInPairedAverage / timesFindInAligned[i].size()) << "\t";
        }
        cout << endl;
    }
    
    cout << "Per alignment operations: " << endl;
    cout << "group\tnum\talignment\ttrace\tpairs" << endl;
    for (int i = 0; i < groups.size(); i++) {
        if (timesAlignment[i].size() == 0) {
            continue;
        }
        double alignmentAverage = 0;
        for (double x : timesAlignment[i]) {
            alignmentAverage += x;
        }
        double traceAverage = 0;
        for (double x : timesTrace[i]) {
            traceAverage += x;
        }
        double pairsAverage = 0;
        for (double x : timesPointPairs[i]) {
            pairsAverage += x;
        }
        
        cout << groups[i] << "\t";
        cout << timesAlignment[i].size() << "\t";
        cout << (alignmentAverage / timesAlignment[i].size()) << "\t";
        cout << (traceAverage / timesTrace[i].size()) << "\t";
        cout << (pairsAverage / timesPointPairs[i].size()) << endl;
    }
    cout << "Test was performed on " << readCounter << " reads." << endl << endl;
}

int GetGroup(int x) {
    int result;
    for (result = 0; result < groups.size(); result++)
        if (x <= groups[result])
            return result;
    throw 1;
}

int ExpandReads(Sequence &a, Sequence &b, vector<Match> &matches, DalignWrapper &dw) {
    AlignedPairsSet alignedPairs;

    int longestAlignment = 0;
    int alignmentsCounter = 0;
    for (Match &match : matches) {
        Utility::StartTiming();
        if (alignedPairs.IsAligned(match)) {
            continue;
        }
        timesFindInAligned[GetGroup(b.GetData().size())].push_back(Utility::GetTimerResult());  

        Alignment al;

        // get compute alignment time
        Utility::StartTiming();
        dw.ComputeAlignment(a, b, match, al);
        double computeAlignmentTime = Utility::GetTimerResult();

        vector<pair<int, int>> pairs;

        // get compute trace time
        Utility::StartTiming();
        al.ComputeTrace();
        double computeTraceTime = Utility::GetTimerResult();   

        // get aligned pairs time
        Utility::StartTiming();
        al.GetAlignedPairs(pairs);
        alignedPairs.MarkAsAligned(pairs);
        double getPairsTime = Utility::GetTimerResult();

        // get index of the group
        int groupIndex = GetGroup(al.GetLengthOnB()); 

        // store value
        timesAlignment[groupIndex].push_back(computeAlignmentTime);
        timesTrace[groupIndex].push_back(computeTraceTime);
        timesPointPairs[groupIndex].push_back(getPairsTime);

        if (al.GetLengthOnB() >= b.GetData().length() - 20)
            break;
//            if (readCounter % 100 == 0)
//                cout << readCounter << endl;
    }
}