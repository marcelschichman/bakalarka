#pragma once
#include "common.h"
#include "Output.h"
using namespace std;

class SAMOutput : public Output {
public:
    SAMOutput(const string& _filename, const Sequence& _genome) :
    Output(_filename, _genome) {
    }
protected:
    virtual bool PrintCustomAlignment(const Alignment& al, const Sequence& read, Direction dir);
    virtual bool PrintCustomHeader(const Sequence& _genome);

};
