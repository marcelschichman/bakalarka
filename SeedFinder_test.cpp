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

int SEED_FINDER_TEST(int argc, char** argv) {
    FASTA fasta("genome.fasta");
    FASTQ fastq("fajny_read.fastq");
    Sequence genome, read;
    fasta >> genome;
    fastq >> read;
    
    vector<Match> forwardSeeds;
    vector<Match> backwardSeeds;
    
    SeedFinder_hashmap fsh(20);
    fsh.CreateIndexFromGenome(genome);
    fsh.GetSeedsWithRead(read, forwardSeeds, backwardSeeds);
    cout << "forward old algorithm: " << endl;
    for (Match &m : forwardSeeds) {
        cout << m.genomePos << " " << m.readPos << " " << m.length << endl;
    }
    cout << endl;
    cout << "backward old algorithm: " << endl;
    for (Match &m : backwardSeeds) {
        cout << m.genomePos << " " << m.readPos << " " << m.length << endl;
    }
    cout << endl;
    
    forwardSeeds.clear();
    backwardSeeds.clear();
    
    
    SeedFinder_hashmap_2bit fsh2b(20);
    fsh2b.CreateIndexFromGenome(genome);
    fsh2b.GetSeedsWithRead(read, forwardSeeds, backwardSeeds);
    cout << "forward new algorithm: " << endl;
    for (Match &m : forwardSeeds) {
        cout << m.genomePos << " " << m.readPos << " " << m.length << endl;
    }
    cout << endl;
    cout << "backward new algorithm: " << endl;
    for (Match &m : backwardSeeds) {
        cout << m.genomePos << " " << m.readPos << " " << m.length << endl;
    }
    cout << endl;
}
