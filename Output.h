#pragma once
#include "common.h"
#include <fstream>
#include "DalignWrapper.h"

using namespace std;

class Output {
public:
    Output(const string &_filename, const Sequence &_genome);
    virtual ~Output();
    bool Init();
    bool AddAlignment(const Alignment &al, const Sequence &read);
protected:
    virtual bool PrintCustomHeader(const Sequence &_genome) = 0;
    virtual bool PrintCustomAlignment(const Alignment &al, const Sequence &read) = 0;

    ofstream os;
private:
    bool isInitialized;
    const string &filename;
    const Sequence &genome;
};
