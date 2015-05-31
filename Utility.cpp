#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <set>

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

void HighErrorRateRegions::AddRegion(const Alignment& left, const Alignment& right) {
    if (left.GetPosOnA().second < right.GetPosOnA().first) {
        regions.push_back(Region());
        Region &reg = regions.back();
        
        reg.pos = make_pair(left.GetPosOnA().second, right.GetPosOnA().first);       
        reg.leftAlignments.push_back(left.GetPosOnA());
        reg.rightAlignments.push_back(right.GetPosOnA());
    }
}

void HighErrorRateRegions::GetHERRegions(int numOccurences, vector<Region>& outRegions) {
    bool inHERRegion = false;
    
    cout << "regions: " << regions.size() << endl;
    
    vector<tuple<int, int, int>> startsAndEnds;
    for (int i = 0; i < (int)regions.size(); i++) {
        startsAndEnds.push_back(make_tuple(regions[i].pos.first, i, 0));
        startsAndEnds.push_back(make_tuple(regions[i].pos.second, i, 1));
    }
    sort(startsAndEnds.begin(), startsAndEnds.end());
    
    set<int> inRegions;
    int maxOverlaping = 0;
    
    outRegions.clear();
    
    for (auto x : startsAndEnds) {
        if (get<2>(x) == 0) {
            inRegions.insert(get<1>(x));
        } else {
            inRegions.erase(get<1>(x));
        }
        
        if (inRegions.size() >= numOccurences) {
            if (!inHERRegion) {
                outRegions.push_back(Region());
                for (int y : inRegions) {
                    outRegions.back().leftAlignments.push_back(regions[y].leftAlignments[0]);
                    outRegions.back().rightAlignments.push_back(regions[y].rightAlignments[0]);
                }
            } else {
                outRegions.back().leftAlignments.push_back(regions[get<1>(x)].leftAlignments[0]);
                outRegions.back().rightAlignments.push_back(regions[get<1>(x)].rightAlignments[0]);
            }
            inHERRegion = true;
        } else {
            inHERRegion = false;
        }
        maxOverlaping = max(maxOverlaping, (int)inRegions.size());
    }
    cout << "inRegions.size(): " << inRegions.size() << endl;
    cout << "max overlaping: " << maxOverlaping << endl;
}

void HighErrorRateRegions::OutputHERRegions(int numOccurences) {
    vector<Region> outRegions;
    int coutner = 0;
    GetHERRegions(numOccurences, outRegions);
    int sumInBadRegions = 0;
    for (auto x : outRegions) {
        int n = x.leftAlignments.size();
        int width = 600;
        int height = 20 * (n + 2);
        cv::Mat img(height, width, CV_8UC3);
        img = cv::Scalar(255, 255, 255);
        
        int startOnGenome = 1234567890;
        int endOnGenome = 0;
        for (auto y : x.leftAlignments) {
            startOnGenome = min(startOnGenome, y.first);
        }
        for (auto y : x.rightAlignments) {
            endOnGenome = max(endOnGenome, y.second);
        }
        int lengthOnGenome = endOnGenome - startOnGenome;
        if (lengthOnGenome <= 0) {
            continue;
        }
        
#define PROJECT2(genomePos, row) cv::Point(((double)(genomePos) - startOnGenome) / lengthOnGenome * (width - 20) + 10, (row) * 20 + 10)
        
        cv::line(img, PROJECT2(startOnGenome, 0), PROJECT2(endOnGenome, 0), cv::Scalar(255, 0, 0), 2);
        for (int i = 0; i < n; i++) {
            cv::line(img, PROJECT2(x.leftAlignments[i].first, i + 1), PROJECT2(x.leftAlignments[i].second, i + 1), cv::Scalar(0, 0, 255), 2);
            cv::line(img, PROJECT2(x.rightAlignments[i].first, i + 1), PROJECT2(x.rightAlignments[i].second, i + 1), cv::Scalar(0, 0, 255), 2);
        }
        
        char buffer[100];
        sprintf(buffer, "bad_regions/%06d.png", coutner++);
        
        // ruller
        cv::Point rullerLeft = PROJECT2(startOnGenome, n + 1);
        cv::Point rullerRight = PROJECT2(startOnGenome + 1000, n + 1);
        cv::Point small(0, 5);
        cv::line(img, rullerLeft, rullerRight, cv::Scalar(0, 0, 0));
        cv::line(img, rullerLeft - small, rullerLeft + small, cv::Scalar(0, 0, 0));
        cv::line(img, rullerRight - small, rullerRight + small, cv::Scalar(0, 0, 0));
        cv::putText(img, "1000", rullerLeft - cv::Point(0, 2), cv::FONT_HERSHEY_PLAIN, 1, cv::Scalar(0, 0, 0));
        
        cv::imwrite(buffer, img);
        sumInBadRegions += n;
    }
    cout << "sum in bad regions: " << sumInBadRegions << endl;
}
