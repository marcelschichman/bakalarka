#include "BlasrOutput.h"

bool BlasrOutput::PrintCustomHeader(const Sequence& _genome) {
    genomeName = _genome.GetId();
    //os << "qName\ttName\tqStrand\ttStrand\tscore\tpercentSimilarity\ttStart\ttEnd\ttLength\tqStart\tqEnd\tqLength\tnCells" << endl;
    os << "qName\ttName\tqStrand\ttStrand\tqLength" << endl;
    
    return true;
}

bool BlasrOutput::PrintCustomAlignment(const Alignment& al, const Sequence& read, Direction dir) {
    os << read.GetId() << "/0_" << read.GetData().length() << "\t";
    os << genomeName << "\t";
    os << "0\t";
    os << (dir == EDIR_FORWARD ? 0 : 1) << "\t";
    os << al.GetLengthOnB() << endl;
}