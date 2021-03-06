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
    bool AddAlignment(const Alignment &al, const Sequence &read, Direction dir);
protected:
    virtual bool PrintCustomHeader(const Sequence &_genome) = 0;
    virtual bool PrintCustomAlignment(const Alignment &al, const Sequence &read, Direction dir) = 0;

    ofstream os;
    const Sequence &genome;
private:
    bool isInitialized;
    const string &filename;
};
