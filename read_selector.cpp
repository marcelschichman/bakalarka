#include "common.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
using namespace std;

ifstream is;
ofstream os;

bool ReadRead(string &a, string &b, string &c, string &d) {
    if (!getline(is, a))
        return false;
    if (!getline(is, b))
        return false;
    if (!getline(is, c))
        return false;
    if (!getline(is, d))
        return false;
    return true;
}

bool WriteRead(const string &a, const string &b, const string &c, const string &d) {
    os << a << endl << b << endl << c << endl << d << endl;
}

int READ_SELECTOR(int argc, char** argv) {
    if (argc < 3) {
        cout << "usage: read_selector inputfile outputfile [number] [min length] [starting postion] [max length]" << endl;
        cout << "negative number means no limit" << endl;
        return 0;
    }
    is.open(argv[1]);
    os.open(argv[2]);
    
    int n = -1, n_orig = -1;
    unsigned int minLength = 0;
    unsigned int maxLength = ~0;
    unsigned int startPos = 0;
    
    if (argc >= 4) {
        n_orig = n = atoi(argv[3]);
    }
    if (argc >= 5) {
        minLength = atoi(argv[4]);
    }
    if (argc >= 6) {
        startPos = atoi(argv[5]);
    }
    if (argc >= 7) {
        maxLength = atoi(argv[6]);
    }
    
    string a, b, c, d;
    while (startPos > 0 && ReadRead(a, b, c, d)) {
        startPos--;
    }
    
    while (n != 0 && ReadRead(a, b, c, d)) {
        if (b.length() >= minLength && b.length() <= maxLength) {
            ostringstream name;
            name << "@READ" << setw(5) << setfill('0') << (n_orig - n);
            WriteRead(name.str(), b, c, d);
            n--;
        }
    }
    cout << "successfully copied " << (n_orig - n) << " reads" << endl;
}
