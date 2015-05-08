#include "common.h"
#include "Sequence.h"
#include "Utility.h"
#include "DalignWrapper.h"

using namespace std;
#define GENOME_POS 1330218
#define READ_POS 571
#define REPETITIONS 50


// 1330218 571

int SIMPLE_TEST(int argc, char** argv) {
    Sequence a("ACTGACGGCTATTACAG"), b("CTGGGCTACATTA");

    DalignWrapper dw;
    dw.SetAligningParameters(0.7, 5, {0.25, 0.25, 0.25, 0.25});
    Alignment al;
    dw.ComputeAlignment(a, b, make_pair(8, 5), al);
    
    al.ComputeTrace();
    vector<pair<int, int>> pairs;
    al.GetAlignedPairs(pairs);
    
    for (auto &p : pairs) {
        cout << p.first << " " << p.second << " " << a.GetData()[p.first] << " " << b.GetData()[p.second] << endl;
    }
    string CIGAR;
    al.GetCigarString(CIGAR);
    cout << CIGAR << endl;
}
