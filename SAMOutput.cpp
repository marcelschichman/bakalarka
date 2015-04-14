#include "SAMOutput.h"

SAMOutput::SAMOutput(const string& filename)
: os(filename.c_str()) {
    os << "@HLAVICKA" << endl;
}

SAMOutput::~SAMOutput() {
}

void SAMOutput::AddAlignment(const Sequence& seq, const Alignment& al) {
    string cigar;
    al.GetCigarString(cigar);
    os << seq.id << " " << seq.data << " " << cigar << endl;
}
