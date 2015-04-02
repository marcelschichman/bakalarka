#pragma once
#include "common.h"
#include <fstream>

using namespace std;


class Sequence {
public:
    Sequence();
    Sequence(const string& _data, string _id = string());
    Sequence(const Sequence& orig);
    virtual ~Sequence();
    
    string data;
    string id;
    
    //TODO: konverzie do roznych formatov
    char* ToDalignFromat();
    
private:
    char* dalignFromat;
};

class FASTA {
public:
    FASTA(const string& _filename);
    
    FASTA& operator>>(Sequence& seq);
private:
    string filename;
};

class FASTQ {
public:
    FASTQ(const string& filename, bool _doReverseRemap = true);
    
    FASTQ& operator>>(Sequence& seq);
    operator bool() const { return isOk; };
private:
    ifstream is;
    bool doReverseRemap;
    bool isOk;
};