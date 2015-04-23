#include "Output.h"

Output::Output(const string& _filename, const Sequence& _genome)
: filename(_filename), genome(_genome) {
}


Output::~Output() {
}

bool Output::Init() {
    os.open(filename.c_str());
    if (os.is_open()) {
        PrintCustomHeader(genome);
        isInitialized = true;
        return true;
    } 
    return false;
}

bool Output::AddAlignment(const Alignment& al, const Sequence& read, Direction dir) {
    if (isInitialized) {
        return PrintCustomAlignment(al, read, dir);
    }
    return false;
}