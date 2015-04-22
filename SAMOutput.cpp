#include "SAMOutput.h"

bool SAMOutput::PrintCustomHeader(const Sequence& _genome) {
    os << "@HLAVICKA" << endl;
}


bool SAMOutput::PrintCustomAlignment(const Alignment& al, const Sequence& read) {
    string cigar;
    al.GetCigarString(cigar);
    os << read.GetId() << " " << read.GetData() << " " << cigar << endl;
}