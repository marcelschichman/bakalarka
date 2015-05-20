#include "BlasrOutput.h"

bool BlasrOutput::PrintCustomHeader(const Sequence& _genome) {
    genomeName = _genome.GetId();
    os << "qName tName qStrand tStrand score percentSimilarity tStart tEnd tLength qStart qEnd qLength nCells" << endl;
    
    return true;
}

bool BlasrOutput::PrintCustomAlignment(const Alignment& al, const Sequence& read, Direction dir) {
    os << read.GetId() << "/0_" << read.GetData().length() << " ";
    os << genomeName << " ";
    os << "0 ";
    os << (dir == EDIR_FORWARD ? 0 : 1) << " ";
    os << "0 ";
    os << "0 ";
    os << al.GetPosOnA().first << " " << al.GetPosOnA().second << " " << genome.GetData().length() << " ";
    os << al.GetPosOnB().first << " " << al.GetPosOnB().second << " " << read.GetData().length() << " ";
    os << 0 << endl;
}