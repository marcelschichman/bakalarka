#pragma once
#include "common.h"
#include "Sequence.h"
#include "SeedFinder.h"
#include "DalignWrapper.h"
#include <chrono>
#include <opencv2/core.hpp>
#include <unordered_set>

typedef vector<pair<vector<cv::Point>, cv::Scalar>> VisualizationData;

class Utility {
////////////////////////////////////////////////////////////////////////////////
// Common
public:
    static int Diagonal(pair<int, int> seed);
    static int AntiDiagonal(pair<int, int> seed);
    
////////////////////////////////////////////////////////////////////////////////
// Match filter
public:
    static void Filter(vector<Match>& matches, int distThreshold, double percentage);
    static void Filter(vector<Match>& forwardMatches, vector<Match>& backwardMatches, int distThreshold, double percentage);
private:
    static void Sort(vector<Match>& matches, int distThreshold);
    
////////////////////////////////////////////////////////////////////////////////
// Timer
public:
    static void StartTiming();
    static double GetTimerResult();
    static double VerbalResult(const string& name);
private:
    static vector<chrono::time_point<chrono::system_clock>> tpts;
    
////////////////////////////////////////////////////////////////////////////////
// Alignment connecting
public:
    static bool AreConnectable(const Alignment& left, const Alignment& right);
    static bool IsCoveredBy(const Alignment& covered, const Alignment& covering);
    static bool IsSignificant(const Alignment& al);
    
    
////////////////////////////////////////////////////////////////////////////////
// Other
public:
};

////////////////////////////////////////////////////////////////////////////////
// Is Aligned
class AlignedPairsSet {
    unordered_set<long long> alignedPairsSet;
public:
    void MarkAsAligned(vector<pair<int, int>> &alignedPairsVector);
    bool IsAligned(Match& seed);
    bool IsAligned(pair<int, int> seed);
};

////////////////////////////////////////////////////////////////////////////////
// Alignment Visualization
class Visualization {
    vector<cv::Point> seeds;
    vector<pair<vector<cv::Point>, cv::Scalar>> alignments;
    string caption;
public:
    Visualization(const string& _caption = "");
    void AddAlignment(Alignment &al, cv::Scalar color = cv::Scalar(0, 0, 255));
    void AddSeeds(vector<Match> &seeds);
    void CreateVisualization(const string &filename = "", int thickness = 1, int outputWidth = -1, int outputHeight = -1, int genomeLength = -1, int readLength = -1);
};

////////////////////////////////////////////////////////////////////////////////
// High Error Rate Regions
class HighErrorRateRegions {
    struct Region {
        pair<int, int> pos;
        vector<pair<int, int>> leftAlignments;
        vector<pair<int, int>> rightAlignments;
    };
    vector<Region> regions;
public:
    void AddRegion(const Alignment &left, const Alignment &right);
    void GetHERRegions(int numOccurences, vector<Region> &outRegions);
    void OutputHERRegions(int numOccurences);
};