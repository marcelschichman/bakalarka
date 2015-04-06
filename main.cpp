#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Timer.h"
#include "DalignWrapper.h"
#include <forward_list>
using namespace std;

int THE_APP(int argc, char** argv) {
    //    Alignment al;
    //    FASTA fasta ("genome.fasta");
    //    Sequence seq;
    //    fasta >> seq;
    //    
    //    SeedFinder_hashmap sf (20);
    //    sf.CreateIndexFromGenome(seq);
    //    
    //    Timer::startTiming();
    //    FASTQ fastq ("pacbio_10kb.fastq");
    //    int i = 0;
    //    while (fastq >> seq) {
    //        Timer::startTiming();
    //        vector<match> result;
    //        sf.GetSeedsWithRead(seq, result);
    //        cout << i++ << " " << seq.id << " seeds: " << result.size() << " time: " << Timer::getTimerResult() << " length: " << seq.data.length() << endl;
    //    }
    //    cout << "time sum: " << Timer::getTimerResult() << endl;


    // get sequences
    Sequence a, b;
    FASTA fasta("genome.fasta");
    FASTQ fastq("pacbio_10kb.fastq");
    fasta >> a;
    //fastq >> b;

    // prepare wrapper
    DalignWrapper dw;
    dw.SetAligningParameters(0.7, 5,{0.25, 0.25, 0.25, 0.25});

    // find seeds
    SeedFinder_hashmap sf(30);
    sf.CreateIndexFromGenome(a);
        
    int readCounter = 0;
    FOR(i, 15) fastq >> b;
    //while (fastq >> b) {
        /*if (b.data.length() < 1000) {
            continue;
        }*/
        cout << readCounter++ << " ";
        
        vector<match> result;
        sf.GetSeedsWithRead(b, result);
        //cout << "total matches: " << result.size() << endl;

        forward_list<match> matches(result.begin(), result.end());
        
        int alignmentsCounter = 0;
        while (!matches.empty()) {
            int i = 0;
            for (auto x : matches) i++;
            Alignment al;
            Timer::startTiming();
            dw.ComputeAlignment(a, b, matches.front(), al);
            Timer::verbalResult("compute alignment");
            cout << "seeds: " << i << " " << al.GetLengthOnB() << " " << b.data.length() << " " << (matches.front().genomePos - matches.front().readPos) << endl;
            cout << matches.front().genomePos << " " << matches.front().readPos << endl;
            matches.pop_front();
            if (al.GetLengthOnB() == b.data.length()) {
                alignmentsCounter++;
                vector<pair<int, int>> pairs;
                Timer::startTiming();
                al.ComputeTrace();
                Timer::verbalResult("Compute trace");
                //al.PrintAlignment("alignment.txt");
                al.GetAlignedPairs(pairs);
                int maxDiagonal = pairs[0].first - pairs[0].second;
                int minDiagonal = maxDiagonal;
                for (auto p : pairs) {
                    minDiagonal = min(minDiagonal, p.first - p.second);
                    maxDiagonal = max(maxDiagonal, p.first - p.second);
                }
                matches.remove_if([minDiagonal, maxDiagonal](match & m) {
                    int diagonal = m.genomePos - m.readPos;
                    return diagonal >= minDiagonal && diagonal <= maxDiagonal;
                });
            }
        }
        cout << "valid alignments: " << alignmentsCounter << endl;
    //}
}

