#pragma once
#include "common.h"
#include "Output.h"
using namespace std;

class SAMOutput : public Output {
public:
    SAMOutput(const string& _filename, const Sequence& _genome) :
    Output(_filename, _genome) {
    }
    virtual bool PrintCustomAlignment(const Alignment& al, const Sequence& read);
    virtual bool PrintCustomHeader(const Sequence& _genome);

};
