#include "common.h"
#include "Sequence.h"
#include "SeedFinder_hashmap.h"
#include "Timer.h"
using namespace std;

int main(int argc, char** argv) {
    FASTA fasta ("genome.fasta");
    Sequence seq;
    fasta >> seq;
    
    SeedFinder_hashmap sf (20);
    sf.CreateIndexFromGenome(seq);
    
    Timer::startTiming();
    FASTQ fastq ("pacbio_10kb.fastq");
    int i = 0;
    while (fastq >> seq) {
        Timer::startTiming();
        vector<match> result;
        sf.GetSeedsWithRead(seq, result);
        cout << i++ << " " << seq.id << " seeds: " << result.size() << " time: " << Timer::getTimerResult() << " length: " << seq.data.length() << endl;
    }
    cout << "time sum: " << Timer::getTimerResult() << endl;
}

