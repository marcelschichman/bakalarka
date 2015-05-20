#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "Utility.h"

vector<chrono::time_point<chrono::system_clock>> Utility::tpts;

int Utility::Diagonal(pair<int, int> seed) {
    return seed.first - (int)(seed.second * 0.92);
}

int Utility::AntiDiagonal(pair<int, int> seed) {
    return seed.first + seed.second;
}

void Utility::Sort(vector<Match>& matches, int distThreshold) {
    vector<int> positions;
    for (Match &match : matches) {
        positions.push_back(Diagonal(match));
    }
    sort(positions.begin(), positions.end());
    
    for (Match &match : matches) {
        match.score = upper_bound(positions.begin(), positions.end(), Diagonal(match) + distThreshold)
                - lower_bound(positions.begin(), positions.end(), Diagonal(match) - distThreshold);
    }
    sort(matches.begin(), matches.end(), [](Match l, Match r) {
        return l.score > r.score;
    });
}

void Utility::Filter(vector<Match>& matches, int distThreshold, double percentage) {
    if (matches.empty()) {
        return;
    }
    Sort(matches, distThreshold);
    
    int highestScore = matches[0].score;
    
    vector<Match> filteredMatches;
    for (Match &match : matches) {
        if (match.score >= highestScore * percentage) {
            filteredMatches.push_back(match);
        }
    }
    matches = filteredMatches;
}

void Utility::Filter(vector<Match>& forwardMatches, vector<Match>& backwardMatches, int distThreshold, double percentage) {
    if (forwardMatches.empty() && backwardMatches.empty()) {
        return;
    }
    Sort(forwardMatches, distThreshold);
    Sort(backwardMatches, distThreshold);

    int highestScore = 0; //max(forwardMatches[0].score, backwardMatches[0].score);
    if (!forwardMatches.empty()) {
        highestScore = max(highestScore, forwardMatches[0].score);
    }
    if (!backwardMatches.empty()) {
        highestScore = max(highestScore, backwardMatches[0].score);
    }
    
    Match ref = {0, 0, 0, 1};
    auto eraseBegin = lower_bound(forwardMatches.begin(), forwardMatches.end(), ref, [](Match l, Match r) {
        return l.score > r.score;
    });
    forwardMatches.erase(eraseBegin, forwardMatches.end());
    
    eraseBegin = lower_bound(backwardMatches.begin(), backwardMatches.end(), ref, [](Match l, Match r) {
        return l.score > r.score;
    });
    backwardMatches.erase(eraseBegin, backwardMatches.end());
}


void Utility::StartTiming() {
    tpts.push_back(chrono::system_clock::now());
}

double Utility::GetTimerResult() {
    double result = ((chrono::duration<double>)(chrono::system_clock::now() - tpts.back())).count();
    tpts.pop_back();
    return result;
}

double Utility::VerbalResult(const string& name) {
    double result = GetTimerResult();
    cout << name << " took " << result << "s" << endl;
    return result;
}

bool Utility::AreConnectable(const Alignment& left, const Alignment& right) {
    pair<int, int> start1, end1, start2, end2;
    start1 = left.GetStartPos();
    end1 = left.GetEndPos();
    start2 = right.GetStartPos();
    end2 = right.GetEndPos();
    
    if (abs(Diagonal(end1) - Diagonal(start2)) < 300 && AntiDiagonal(start1) < AntiDiagonal(start2) &&
        AntiDiagonal(end1) < AntiDiagonal(end2) && AntiDiagonal(end1) - AntiDiagonal(start2) < 500) {
        return true;
    }
    return false;
}

bool Utility::IsCoveredBy(const Alignment& covered, const Alignment& covering) {
    pair<int, int> start1, end1, start2, end2;
    start1 = covered.GetStartPos();
    end1 = covered.GetEndPos();
    start2 = covering.GetStartPos();
    end2 = covering.GetEndPos();
    if ((Diagonal(start1) == Diagonal(start2) || Diagonal(start1) == Diagonal(end2)) &&
        (Diagonal(end1) == Diagonal(start2) || Diagonal(end1) == Diagonal(end2)) &&
        AntiDiagonal(start1) >= AntiDiagonal(start2) && AntiDiagonal(end1) <= AntiDiagonal(end2)) {
        return true;
    }
    return false;
}

bool Utility::IsSignificant(const Alignment& al) {
    return al.GetLengthOnB() >= 200;
}

void AlignedPairsSet::MarkAsAligned(vector<pair<int, int> >& alignedPairsVector) {
    long long previous = -1;
    for (auto &p : alignedPairsVector) {
        long long numValue = (((long long)p.first >> 3) << 32) + (p.second >> 3);
        if (numValue != previous) {
            alignedPairsSet.insert(numValue);
            previous = numValue;
        }
    }
}

Visualization::Visualization(const string& _caption)
: caption(_caption){

}


bool AlignedPairsSet::IsAligned(Match& seed) {
    long long numValue = (((long long)(seed.genomePos + seed.length / 2) >> 3) << 32) + ((seed.readPos + seed.length / 2) >> 3);
    return alignedPairsSet.find(numValue) != alignedPairsSet.end();
}

bool AlignedPairsSet::IsAligned(pair<int,int> seed) {
    long long numValue = (((long long)(seed.first) >> 3) << 32) + ((seed.second) >> 3);
    return alignedPairsSet.find(numValue) != alignedPairsSet.end();
}

void Visualization::CreateVisualization(const string &filename, int thickness, int outputWidth, int outputHeight, int genomeLength, int readLength) {
    int minGenome = 1234567890;
    int maxGenome = 0;
    int minRead = 1234567890;
    int maxRead = 0;
    if (genomeLength == -1) {
        for (auto al : alignments) {
            for (auto p : al.first) {
                minGenome = min(minGenome, p.x);
                maxGenome = max(maxGenome, p.x);
            }
        }
        for (auto p : seeds) {
            minGenome = min(minGenome, p.x);
            maxGenome = max(maxGenome, p.x);
        }
    } else {
        minGenome = 0;
        maxGenome = genomeLength - 1;
    }
    
    if (readLength == -1) {
        for (auto al : alignments) {
            for (auto p : al.first) {
                minRead = min(minRead, p.y);
                maxRead = max(maxRead, p.y);
            }
        }
        for (auto p : seeds) {
            minRead = min(minRead, p.y);
            maxRead = max(maxRead, p.y);
        }
    } else {
        minRead = 0;
        maxRead = readLength - 1;
    }
    
    if (maxRead - minRead <= 0 || maxGenome - minGenome <= 0) {
        return;
    }
    
    if (outputWidth < 0) {
        outputWidth = maxGenome - minGenome;
    }
    if (outputHeight < 0) {
        outputHeight = maxRead - minRead;
    }
    
    cv::Mat img(outputHeight, outputWidth, CV_8UC3);
    img = cv::Scalar(255, 255, 255);
    if (caption != "") {
        cv::putText(img, caption, cv::Point(0, 15), CV_FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 0, 0));
    }
    
    cv::Point topLeft(minGenome, minRead);
    
#define PROJECT(point) cv::Point(((double)(point.x) - topLeft.x) / (maxGenome - minGenome) * (outputWidth), \
    ((double)(point.y) - topLeft.y) / (maxRead - minRead) * (outputHeight))    
    
    for (auto al : alignments) {
        if (al.first.size() == 0) {
            continue;
        }
        
        cv::Point prev = PROJECT(al.first[0]);
        for (auto p : al.first) {
            cv::line(img, prev, PROJECT(p), al.second, thickness);
            prev = PROJECT(p);
        }
    }
    
    for (auto p : seeds) {
        cv::circle(img, PROJECT(p), thickness, cv::Scalar(0, 0, 0), -1);
    }
    
    if (filename == "") {
        cv::imshow("Alignments", img);
        cv::waitKey(0);
    } else {
        cv::imwrite(filename, img);
    }
}

void Visualization::AddAlignment(Alignment &al, cv::Scalar color) {
    pair<vector<cv::Point>, cv::Scalar> alv;
    alv.second = color;
    vector<pair<int, int>> pairs;
    al.GetAlignedPairs(pairs);
    for (auto p : pairs) {
        alv.first.push_back(cv::Point(p.first, p.second));
    }
    alignments.push_back(alv);
}

void Visualization::AddSeeds(vector<Match>& seeds) {
    for (auto seed : seeds) {
        this->seeds.push_back(cv::Point(seed.genomePos, seed.readPos));
    }
}
