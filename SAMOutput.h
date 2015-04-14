#pragma once
#include "common.h"
#include <fstream>
#include "DalignWrapper.h"
using namespace std;

class SAMOutput {
public:
    SAMOutput(const string &filename);
    void AddAlignment(const Sequence& seq, const Alignment& al);
    virtual ~SAMOutput();
private:
    ofstream os; 
};
