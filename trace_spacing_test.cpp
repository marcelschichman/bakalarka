#include "common.h"
#include "Sequence.h"
#include "Timer.h"
#include "DalignWrapper.h"

using namespace std;
#define GENOME_POS 1330218
#define READ_POS 571
#define REPETITIONS 50

Sequence a, b;

void test(int spacing) {
    DalignWrapper dw;
    Alignment al;
    dw.SetAligningParameters(0.7, spacing,{0.25, 0.25, 0.25, 0.25});

    double compAlTimeSum = 0;
    double compTrTimeSum = 0;

    FOR(i, REPETITIONS) {
        Timer::startTiming();
        dw.ComputeAlignment(a, b, make_pair(GENOME_POS, READ_POS), al);
        compAlTimeSum += Timer::getTimerResult();
        ;

        Timer::startTiming();
        al.ComputeTrace();
        compTrTimeSum += Timer::getTimerResult();
        ;
    }
    cout << spacing << "\t" << compAlTimeSum / REPETITIONS << "\t" << compTrTimeSum / REPETITIONS << endl;
}

// 1330218 571

int TRACE_SPACING_TEST(int argc, char** argv) {
    cout << "*** Trace Spacing test ***" << endl;
    cout << "The purpose of this test is to decide whether execution time of Local_alignment and Compute_trace_MID functions are influenced by trace spacing parameter." << endl;

    FASTA fasta("genome.fasta");
    FASTQ fastq("pacbio_10kb.fastq");
    fasta >> a;
    FOR(i, 15) fastq >> b;

    cout << "spacing\talignment\ttrace\n";
    test(1);
    test(2);
    test(5);
    test(10);
    test(20);
    test(50);
    test(100);
    test(200);
    test(500);
    test(1000);
}
