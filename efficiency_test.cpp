#include "common.h"
#include "Sequence.h"
#include <map>
#include <sstream>
using namespace std;

struct posInfo {
    int start, end;
    Direction dir;
    double percentageFound;
    string genomeRow, readRow;
};

#define ACCEPT_DIST 300

int EFFICIENCY_TEST(int argc, char** argv) {
    string refFilename = "sd_0001.maf";
    string blsFilename = "pbsim_out.txt";
    
    string dummy, name;
    map<string, posInfo> positions;
    
    if (argc == 3) {
        refFilename = argv[1];
        blsFilename = argv[2];
    }
    
    ifstream is(refFilename);
    while (is.is_open() && !is.eof()) {
        string genomeRow, readRow;
        int startOnGenome, lengthOnGenome;
        int genomeLength;
        Direction dir;
        
        getline(is, dummy);
        getline(is, genomeRow);
        getline(is, readRow);
        getline(is, dummy);
        
        stringstream ss(genomeRow);
        
        ss >> dummy >> dummy >> startOnGenome >> lengthOnGenome >> dummy >> genomeLength >> dummy;
        
        string plusMinus;
        string readName;
        stringstream ss2(readRow);
        ss2 >> dummy >> readName >> dummy >> dummy >> plusMinus >> dummy >> dummy;
        
        if (plusMinus == "-") {
            startOnGenome = genomeLength - (startOnGenome + lengthOnGenome);
            
        }
        
        positions[readName] = {startOnGenome, startOnGenome + lengthOnGenome, plusMinus == "+" ? EDIR_FORWARD : EDIR_BACKWARD, 0, genomeRow, readRow};
    }
    
    ifstream is2(blsFilename);
    getline(is2, dummy);
    
    while (is2.is_open() && !is2.eof()) {
        string readName;
        int orientation;
        int startOnGenome, endOnGenome;
        int startOnRead, endOnRead, readLength;
        
        is2 >> readName >> dummy >> dummy >> orientation >> dummy >> dummy >> startOnGenome >> endOnGenome;
        is2 >> dummy >> startOnRead >> endOnRead >> readLength >> dummy;
        
        readName = readName.substr(0, readName.find("/"));
        
        if (positions.find(readName) == positions.end()) {
            continue;
        }
        
        posInfo &pos = positions[readName];
        
        startOnGenome -= startOnRead;
        endOnGenome += (readLength - endOnRead);
        Direction dir = orientation ? EDIR_BACKWARD : EDIR_FORWARD;
        
        if (abs(pos.start - startOnGenome) < ACCEPT_DIST && abs(pos.end - endOnGenome) < ACCEPT_DIST && pos.dir == dir) {
            pos.percentageFound = max(pos.percentageFound, (double)(endOnRead - startOnRead) / readLength);
        }
    }
    
    // results
    cout << "reads: " << positions.size() << endl;
    vector<double> thresholds = {1, 0.95, 0.90, 0.8, 0.7, 0.5, 0.3, 0.1};
    
    for (double threshold : thresholds) {
        cout << "aligned (" << (threshold * 100) << "%): " << count_if(positions.begin(), positions.end(), [threshold](pair<string, posInfo> p) {
            return p.second.percentageFound >= threshold;
        }) << endl;
    }
    
    ofstream os("not_aligned.maf");
    for (auto pi : positions) {
        if (pi.second.percentageFound < 0.1) {
            os << "a" << endl << pi.second.genomeRow << endl << pi.second.readRow << endl << endl;
        }
    }
    os.flush();
}
