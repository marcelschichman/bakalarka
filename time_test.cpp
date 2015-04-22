#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Timer.h"
#include "DalignWrapper.h"
#include "SAMOutput.h"
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

#define NUM_READS 100000

int TIME_TEST(int argc, char** argv) {
    //header
    cout << "*** Execution time test ***" << endl;
    cout << "This tests tests how much the length of the alignment influences the execution time of function calls." << endl;
    cout << "One time operations: " << endl;
    
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
    
    Timer::startTiming();
    fasta >> a;
    Timer::verbalResult("Reading genome");
    cout << "Genome length: " << a.GetData().length() << endl;
    
    //SAMOutput soutput(outputFilename);

    // prepare wrapper
    Timer::startTiming();
    DalignWrapper dw;
    dw.SetAligningParameters(0.70, 5, {0.25, 0.25, 0.25, 0.25});
    Timer::verbalResult("Preparing DALIGN");

    // find seeds
    Timer::startTiming();
    SeedFinder_hashmap sf(20);
    sf.CreateIndexFromGenome(a);
    Timer::verbalResult("Preparing seed finder");
        
    int readCounter = 0;
    int progressBarPos = 0;
    
    while (fastq >> b) {
        // get seed finder time
        Timer::startTiming();
        vector<Match> matches;
        sf.GetSeedsWithRead(b, matches);
        timesFindSeeds[GetGroup(b.GetData().size())].push_back(Timer::getTimerResult());        
        
        set<pair<int, int>> alignedPairs;
        
        int longestAlignment = 0;
        int alignmentsCounter = 0;
        for (Match &match : matches) {
            Timer::startTiming();
            if (alignedPairs.find(make_pair(match.genomePos+8, match.readPos+8)) != alignedPairs.end()) {
                continue;
            }
            timesFindInAligned[GetGroup(b.GetData().size())].push_back(Timer::getTimerResult());  
            
            Alignment al;
            
            // get compute alignment time
            Timer::startTiming();
            dw.ComputeAlignment(a, b, match, al);
            double computeAlignmentTime = Timer::getTimerResult();
            
            vector<pair<int, int>> pairs;
            
            // get compute trace time
            Timer::startTiming();
            al.ComputeTrace();
            double computeTraceTime = Timer::getTimerResult();   
            
            // get aligned pairs time
            Timer::startTiming();
            al.GetAlignedPairs(pairs);
            alignedPairs.insert(pairs.begin(), pairs.end());
            double getPairsTime = Timer::getTimerResult();
            
            // get index of the group
            int groupIndex = GetGroup(al.GetLengthOnB()); 
            
            // store value
            timesAlignment[groupIndex].push_back(computeAlignmentTime);
            timesTrace[groupIndex].push_back(computeTraceTime);
            timesPointPairs[groupIndex].push_back(getPairsTime);
            
//            if (readCounter % 100 == 0)
//                cout << readCounter << endl;
        }
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