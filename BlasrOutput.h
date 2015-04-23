#pragma once
#include "Output.h"

class BlasrOutput : public Output  {
public:
    BlasrOutput(const string& _filename, const Sequence& _genome) :
    Output(_filename, _genome) {
    }
protected:
    virtual bool PrintCustomAlignment(const Alignment& al, const Sequence& read, Direction dir);
    virtual bool PrintCustomHeader(const Sequence& _genome);
    
    string genomeName;
};
