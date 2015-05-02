#include "MatchFilter.h"
#include <map>
#include <algorithm>

MatchFilter::MatchFilter(vector<Match>& _matches)
: matches(_matches) {
    
}

void MatchFilter::Sort(int nearDist) {
    vector<int> positions;
    for (Match &match : matches) {
        positions.push_back(Diagonal(match));
    }
    sort(positions.begin(), positions.end());
    
    for (Match &match : matches) {
        match.score = upper_bound(positions.begin(), positions.end(), Diagonal(match) + nearDist)
                - lower_bound(positions.begin(), positions.end(), Diagonal(match) - nearDist);
    }
    sort(matches.begin(), matches.end(), [](Match l, Match r) {
        return l.score > r.score;
    });
}

int MatchFilter::Diagonal(Match& match) {
    return match.genomePos - match.readPos;
}

void MatchFilter::Filter(double remove) {
    if (matches.empty()) {
        return;
    }
    int highestScore = matches[0].score;
    
    vector<Match> filteredMatches;
    for (Match &match : matches) {
        if (match.score >= highestScore * remove) {
            filteredMatches.push_back(match);
        }
    }
    matches = filteredMatches;
}
