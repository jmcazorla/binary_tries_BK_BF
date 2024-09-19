#include <iostream>
#include <vector>
#include <stack>
#include <chrono>
#include <thread>
#include <mutex>
#include <queue>
#include <tuple>
#include <math.h>
#include "parallel_for.hpp"

using namespace std;


template<class trieType>
void runsAND(vector<trieType> &Ts, uint64_t nTries, uint64_t &maxLevel,
        uint64_t currLevel, uint64_t roots[], bool activeTries[], 
        uint64_t prefix, vector<uint64_t> &r, uint64_t &delta){
    if (currLevel == maxLevel) {
        delta++;
        r.push_back(prefix);
        return;
    }
    uint64_t i;
    bool tempActiveTries[16];
    uint64_t andResult = 0b11;
    uint64_t orResult = 0b00;
    for (i = 0; i < nTries; ++i) {
        if (activeTries[i]) {
            uint64_t node_i = (currLevel == maxLevel -1) ?
                                Ts[i].getNode2(roots[i]):
                                Ts[i].getNode1(roots[i]);
            if (node_i) {
                tempActiveTries[i] = true;
                andResult &= node_i;
            } else {
                tempActiveTries[i] = false;
            }
            orResult |= node_i;
        } else {
            tempActiveTries[i] = false;
        }
    }

    if (andResult == 0b00 && orResult != 0b00) {
        delta += 2;
        return;
    }

    if (orResult = 0b00) {
        delta++;
        uint64_t below = prefix;
        uint64_t range = ((uint64_t)1 << (maxLevel - currLevel))- 1;
        uint64_t above = prefix | range;
        for (uint64_t i = below; i <= above; ++i) {
            r.push_back(i);
        }
        return;
    }

    bool existLeft, existRight;
    uint64_t nextLevel = currLevel + 1;
    if(andResult == 0b10) {
        delta++;
        uint64_t leftNodes[16];
        uint64_t leftResult = prefix;
        for (i = 0; i < nTries; ++i) {
            if (tempActiveTries[i] && currLevel != maxLevel - 1) 
                leftNodes[i] = Ts[i].getLeftChild(roots[i]);
        }
        runsAND(Ts, nTries, maxLevel, nextLevel, leftNodes,
                        tempActiveTries, leftResult, r, delta);
    }
    else if (andResult == 0b01) {
        delta++;
        uint64_t rightNodes[16];
        uint64_t rightResult = (prefix | (1ULL << (maxLevel- currLevel - 1)));
        for (i = 0; i < nTries; ++i) {
            if (tempActiveTries[i] && currLevel != maxLevel -1)
                rightNodes[i] = Ts[i].getRightChild(roots[i]);
        }
        runsAND(Ts, nTries, maxLevel, nextLevel, rightNodes,
                        tempActiveTries, rightResult, r, delta);

    }
    else if (andResult == 0b11) {
        uint64_t leftNodes[16];
        uint64_t rightNodes[16];
        uint64_t leftResult = prefix;
        uint64_t rightResult = (prefix | (1ULL << (maxLevel- currLevel - 1)));;
        for (i = 0; i < nTries; ++i) {
            if (tempActiveTries[i] && currLevel != maxLevel - 1) {
                uint64_t leftNode = Ts[i].getLeftChild(roots[i]);
                leftNodes[i] = leftNode;
                rightNodes[i] = leftNode + 1;
            }
        }
        runsAND(Ts, nTries, maxLevel, nextLevel, leftNodes,
                        tempActiveTries, leftResult, r, delta);
        runsAND(Ts, nTries, maxLevel, nextLevel, rightNodes,
                        tempActiveTries, rightResult, r, delta);
    }
    return;
};


template<class trieType>
bool AND(vector<trieType> &Ts, uint64_t nTries, uint64_t &maxLevel,
        uint64_t currLevel, uint64_t roots[], uint64_t prefix, vector<uint64_t> &r, uint64_t &delta){
    if (currLevel == maxLevel) {
        r.push_back(prefix);
        delta++;
        return true;
    }
    uint64_t i;
    uint64_t andResult = 0b11;
    for (i = 0; i < nTries; ++i) {
        uint64_t node_i = (currLevel == maxLevel -1) ?
                            Ts[i].getNode2(roots[i]):
                            Ts[i].getNode1(roots[i]);
        andResult &= node_i;
        
    }
    // Can't go any further down in that branch.
    if (andResult == 0b00) {
        delta += 2;
        return false;
    }

    bool existLeft, existRight;
    uint64_t nextLevel = currLevel + 1;
    if(andResult == 0b10) {
        delta++;
        uint64_t leftNodes[16];
        uint64_t leftResult = prefix;
        for (i = 0; i < nTries; ++i) {
            if (currLevel != maxLevel - 1) 
                leftNodes[i] = Ts[i].getLeftChild(roots[i]);
        }
        existLeft = AND(Ts, nTries, maxLevel, nextLevel, leftNodes, leftResult, r, delta);
    }
    else if (andResult == 0b01) {
        delta++;
        uint64_t rightNodes[16];
        uint64_t rightResult = (prefix | (1ULL << (maxLevel- currLevel - 1)));
        for (i = 0; i < nTries; ++i) {
            if (currLevel != maxLevel -1)
                rightNodes[i] = Ts[i].getRightChild(roots[i]);
        }
        existRight = AND(Ts, nTries, maxLevel, nextLevel, rightNodes, rightResult, r, delta);

    }
    else if (andResult == 0b11) {
        uint64_t leftNodes[16];
        uint64_t rightNodes[16];
        uint64_t leftResult = prefix;
        uint64_t rightResult = (prefix | (1ULL << (maxLevel- currLevel - 1)));
        for (i = 0; i < nTries; ++i) {
            if (currLevel != maxLevel - 1) {
                uint64_t leftNode = Ts[i].getLeftChild(roots[i]);
                leftNodes[i] = leftNode;
                rightNodes[i] = leftNode + 1;
            }
        }
        existLeft = AND(Ts, nTries, maxLevel, nextLevel, leftNodes, leftResult, r, delta);
        existRight = AND(Ts, nTries, maxLevel, nextLevel, rightNodes, rightResult, r, delta);
    }
    if (existLeft || existRight)
        return true;
    return false;
};


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template <class trieType>
void DeltaTrieCertificate(vector<trieType> &Bs, vector<uint64_t> &result, uint64_t &delta, bool runs){
    uint64_t max_level = Bs[0].getHeight();
    uint64_t n_tries = Bs.size();
    // max 16 relations
    bool activeTries[16] = { true, true, true, true,
                             true, true, true, true,
                             true, true, true, true,
                             true, true, true, true };
    uint64_t roots[16] = { 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0 };
    
    uint64_t curr_level = 0;
   
    result.reserve(1000000);
 

    if (runs)
        runsAND(Bs, n_tries, max_level, 0, roots, activeTries, 0, result, delta);
    else 
        AND(Bs, n_tries, max_level, 0, roots, 0, result, delta);

    
    return;
}

