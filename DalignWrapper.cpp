#include "DalignWrapper.h"

DalignWrapper::DalignWrapper()
: workData(NULL), alignSpec(NULL) {
}

DalignWrapper::~DalignWrapper() {
    FreeDalignData();
}

void DalignWrapper::SetAligningParameters(float corelation, int traceSpacing, const array<float, 4>& frequencies) {
    FreeDalignData();
    float freq[4] = {frequencies[0], frequencies[1], frequencies[2], frequencies[3]};
    alignSpec = dalign::New_Align_Spec(corelation, traceSpacing, freq);
    workData = dalign::New_Work_Data();
}

bool DalignWrapper::ComputeAlignment(Sequence& A, Sequence& B, pair<int, int> seed, Alignment& al) {
    if (workData == NULL || alignSpec == NULL) {
        return false;
    }
    al.Prepare(A, B, workData, dalign::Trace_Spacing(alignSpec));
    dalign::Local_Alignment(&(al.alignment), workData, alignSpec, seed.first, seed.first, seed.second);
    al.status = Alignment::AS_ALIGNMENT;
    return true;
}

void DalignWrapper::FreeDalignData() {
    if (workData != NULL) {
        dalign::Free_Work_Data(workData);
        workData = NULL;
    }
    if (alignSpec != NULL) {
        dalign::Free_Align_Spec(alignSpec);
        alignSpec = NULL;
    }
}

Alignment::Alignment()
: status(AS_EMPTY) {
}

void Alignment::Prepare(Sequence& _A, Sequence& _B, dalign::Work_Data* _workData, int _traceSpacing) {
    status = AS_PREPARED;
    workData = _workData;
    traceSpacing = _traceSpacing;
    alignment.path = &path;
    alignment.aseq = _A.ToDalignFromat();
    alignment.alen = _A.data.length();
    alignment.bseq = _B.ToDalignFromat();
    alignment.blen = _B.data.length();
    alignment.flags = 0;
}

int Alignment::GetLengthOnA() {
    if (status < AS_ALIGNMENT) {
        return -1;
    }

    return path.aepos - path.abpos;
}

int Alignment::GetLengthOnB() {
    if (status < AS_ALIGNMENT) {
        return -1;
    }

    return path.bepos - path.bbpos;
}

bool Alignment::ComputeTrace() {
    if (status < AS_ALIGNMENT) {
        return false;
    }

    dalign::Compute_Trace_MID(&alignment, workData, traceSpacing);
    status = AS_TRACE;
    return true;
}

bool Alignment::PrintAlignment(const string& filename) {
    if (status < AS_TRACE) {
        return false;
    }

    FILE *alignmentFile = fopen(filename.c_str(), "w");
    Print_Alignment(alignmentFile, &alignment, workData, 10, 80, 5, 1, 10);
    return true;
}

bool Alignment::GetAlignedPairs(vector<pair<int, int> >& pairs) {
    pairs.clear();
    if (status < AS_TRACE) {
        return false;
    }
    int p, c; /* Output columns of alignment til reach trace end */
    int i = path.abpos + 1;
    int j = path.bbpos + 1;
    // copy-pasterino
    for (c = 0; c < path.tlen; c++)
        if ((p = ((int*)path.trace)[c]) < 0) {
            p = -p;
            while (i != p) {
                pairs.push_back(make_pair(i, j));
                i += 1;
                j += 1;
            }
            j += 1;
        } else {
            while (j != p) {
                pairs.push_back(make_pair(i, j));
                i += 1;
                j += 1;
            }
            i += 1;
        }
    p = path.aepos;
    while (i <= p) {
        pairs.push_back(make_pair(i, j));
        i += 1;
        j += 1;
    }
    
    return true;
}