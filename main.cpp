#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Timer.h"
#include "DalignWrapper.h"
using namespace std;

int main(int argc, char** argv) {
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
    
    Alignment al;
    DalignWrapper dw;
    dw.SetAligningParameters(0.7, 2, {0.25, 0.25, 0.25, 0.25});
    
    Sequence a ("AAACCCTTTGGGAAAC");
    Sequence b ("CCCTGGG");
    dw.ComputeAlignment(a, b, make_pair(8, 3), al);
    cout << "length on a: " << al.GetLengthOnA() << endl;
    cout << "length on b: " << al.GetLengthOnB() << endl;
    cout << "aligned bases: " << endl;
    vector<pair<int, int>> pairs;
    al.ComputeTrace();
    al.GetAlignedPairs(pairs);
    for (auto p : pairs) {
        cout << p.first << " " << p.second << endl;
    }
}

