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
#include "intersection_brute_force.hpp"
#include "barbay_and_kenyon.hpp"

using namespace std;

template <class HybridStructure>
void partialAND_HS(vector<HybridStructure> &setHS, uint16_t n_tries, uint64_t max_level, uint64_t curr_level, 
                uint64_t cut_level, uint64_t roots[], bool activeHS[], uint64_t partial_int,
                vector <uint64_t> &r, vector<uint64_t> &partial_ints, vector<uint64_t*> &threads_roots,
                vector<bool*> &threads_activeTries) {
    if (curr_level == cut_level) {
        partial_ints.push_back(partial_int);
        uint64_t* rts = new uint64_t[16];
        bool* atrs;
        atrs = new bool[16];
        for (uint64_t i = 0; i < n_tries; ++i){
            rts[i] = roots[i];
        }
        threads_roots.push_back(rts);        
        threads_activeTries.push_back(atrs);
        return;
    }
    
    bool tempActiveTries [16]; 
    uint64_t result = 0b11;
    uint64_t node00 = 0b00;
    for (uint16_t i = 0;  i < n_tries; ++i) {
        uint64_t node_i = (curr_level == max_level - 1) ? 
                                setHS[i].getBinTrie()->getNode2(roots[i]):
                                setHS[i].getBinTrie()->getNode1(roots[i]);
        result &= node_i;        
    }

    if (result == 0b00) 
        return;
    if (result == 0b10) {
        uint64_t left_nodes[16];
        uint64_t leftResult = partial_int;
        for (uint64_t i = 0; i < n_tries; ++i) {
            left_nodes[i] = setHS[i].getBinTrie()->getLeftChild(roots[i]);
        }
        partialAND_HS(setHS, n_tries, max_level, curr_level + 1, cut_level, 
            left_nodes, tempActiveTries, leftResult, r, partial_ints,
            threads_roots, threads_activeTries);
    }
    else if (result == 0b01){
        uint64_t right_nodes[16];
        uint64_t rightResult = partial_int;
        rightResult = (rightResult | (1ULL << (max_level- curr_level - 1)));
        for (uint64_t i = 0; i < n_tries; ++i) {
    
            right_nodes[i] = setHS[i].getBinTrie()->getRightChild(roots[i]);
        }
        partialAND_HS(setHS, n_tries, max_level, curr_level + 1, cut_level,
            right_nodes, tempActiveTries, rightResult, r, partial_ints,
            threads_roots, threads_activeTries);
    }
    else if (result == 0b11) {
        uint64_t left_nodes[16];
        uint64_t right_nodes[16];
        uint64_t leftResult = partial_int;
        uint64_t rightResult = partial_int;
        rightResult = (rightResult | (1ULL << (max_level- curr_level - 1)));
        for(uint64_t i = 0; i < n_tries; ++i) {
            uint64_t left_node = setHS[i].getBinTrie()->getLeftChild(roots[i]);
            left_nodes[i]  = left_node;
            right_nodes[i] = left_node + 1;        
        }
        // Left Childs
        partialAND_HS(setHS, n_tries, max_level, curr_level + 1, cut_level, 
            left_nodes, tempActiveTries, leftResult, r, partial_ints,
            threads_roots, threads_activeTries);
        // Right Childs
        partialAND_HS(setHS, n_tries, max_level, curr_level + 1, cut_level,
            right_nodes, tempActiveTries, rightResult, r, partial_ints,
            threads_roots, threads_activeTries);
    }
}

template <class HybridStructure, uint8_t low_part_size>
bool AND_High_BTrie_Log(vector<HybridStructure> &setHS, uint64_t nTries, uint64_t &maxLevel,
    uint64_t currLevel, uint64_t roots[], uint64_t prefix, vector<uint64_t> &r){


    uint64_t i;
    uint64_t andResult = 0b11;

    for (i = 0; i < nTries; ++i) {
        uint64_t node_i = (currLevel == maxLevel -1) ?
                            setHS[i].getBinTrie()->getNode2(roots[i]):
                            setHS[i].getBinTrie()->getNode1(roots[i]);
        andResult &= node_i;
    }

    if (andResult == 0b00)
        return false;


    bool existLeft, existRight;
    uint64_t nextLevel = currLevel + 1;
    if(andResult == 0b10) {
        uint64_t leftNodes[16];
        uint64_t leftR_pos[16];
        uint64_t leftResult = prefix;

        if (currLevel != maxLevel - 1) {
            for (i = 0; i < nTries; ++i) {
                leftNodes[i] = setHS[i].getBinTrie()->getLeftChild(roots[i]);
            }
        }
        else {
            uint64_t firstL[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            uint64_t secondL[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            for (i = 0; i < nTries; ++i) {
                uint64_t pos =  setHS[i].getBinTrie()->getChildInLastLevel(roots[i]);

                firstL[i] = pos <= 1 ? 0 : setHS[i].selectFlags(pos -1) + 1;
                secondL[i] = setHS[i].selectFlags(pos);
            }             
            uint64_t min = secondL[0] - firstL[0];
            for (i = 1; i < nTries; ++i) {                
                uint64_t lengthL = secondL[i] - firstL[i];
                if (lengthL < min)
                    min = lengthL;
            }
            if (min > 0)
            {
                r.push_back(min);
            } 
            // vector<uint64_t> intersection_AND_HS; 
            // barbayKenyonHS<HybridStructure, low_part_size>(setHS, firstL, secondL, nTries, leftResult, r);
            // r.insert(r.end(), intersection_AND_HS.begin(), intersection_AND_HS.end());
            return true;
        }
        existLeft = AND_High_BTrie_Log<HybridStructure, low_part_size>(setHS, nTries, maxLevel, nextLevel, leftNodes, leftResult, r);
    }
    else if (andResult == 0b01) {
        uint64_t rightNodes[16];
        uint64_t rightR_pos[16];
        uint64_t rightResult = (prefix | (1ULL << (maxLevel- currLevel - 1)));
        if (currLevel != maxLevel -1) {
            for (i = 0; i < nTries; ++i) {
                rightNodes[i] = setHS[i].getBinTrie()->getRightChild(roots[i]);
            }
        }
        else {
            uint64_t firstR[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            uint64_t secondR[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            for (i = 0; i < nTries; ++i) {
                uint64_t pos =  setHS[i].getBinTrie()->getChildInLastLevel(roots[i]);

                firstR[i] = pos <= 1 ? 0 : setHS[i].selectFlags(pos) + 1;
                secondR[i] = setHS[i].selectFlags(pos+1);

            }
            uint64_t min = secondR[0] - firstR[0];
            for (i = 1; i < nTries; ++i) {                
                uint64_t lengthR = secondR[i] - firstR[i];
                if (lengthR < min)
                    min = lengthR;
            }
            if (min > 0)
            {
                r.push_back(min);
            } 
            // vector<uint64_t> intersection_AND_HS_R; 
            // barbayKenyonHS<HybridStructure, low_part_size>(setHS, firstR, secondR, nTries, rightResult, r);
            // r.insert(r.end(), intersection_AND_HS_R.begin(), intersection_AND_HS_R.end());
            return true;
        }
        existRight = AND_High_BTrie_Log<HybridStructure, low_part_size>(setHS, nTries, maxLevel, nextLevel, rightNodes, rightResult, r);

    }
    else if (andResult == 0b11) {
        uint64_t leftNodes[16];
        uint64_t rightNodes[16];
        uint64_t rightR_pos[16];
        uint64_t leftR_pos[16];
        uint64_t leftResult = prefix;
        uint64_t rightResult = (prefix | (1ULL << (maxLevel- currLevel - 1)));
        if (currLevel != maxLevel - 1) {
            for (i = 0; i < nTries; ++i) {
                uint64_t leftNode = setHS[i].getBinTrie()->getLeftChild(roots[i]);
                leftNodes[i] = leftNode;
                rightNodes[i] = leftNode + 1;
            }
        }
        else {
            uint64_t firstL[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            uint64_t secondL[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

            uint64_t firstR[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            uint64_t secondR[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            for (i = 0; i < nTries; ++i) {
                uint64_t pos =  setHS[i].getBinTrie()->getChildInLastLevel(roots[i]);

                firstL[i] = pos <= 1 ? 0 : setHS[i].selectFlags(pos -1) + 1;
                secondL[i] = setHS[i].selectFlags(pos);

                firstR[i] = setHS[i].selectFlags(pos) + 1;
                secondR[i] = setHS[i].selectFlags(pos+1); 
            }
            uint64_t minL = secondL[i] - firstL[i];
            uint64_t minR = secondR[i] - firstR[i];
            for (i = 1; i < nTries; ++i) {                
                uint64_t lengthR = secondL[i] - firstL[i];
                if (lengthR < minR)
                    minR = lengthR;
                uint64_t lengthL = secondL[i] - firstL[i];
                if (lengthL < minL)
                    minL = lengthL;
            }
            if (minL > 0)
            {
                r.push_back(minL);
            } 
            if (minL > 0)
            {
                r.push_back(minR);
            }               

            
            // vector<uint64_t> intersection_AND_HS; 
            // barbayKenyonHS<HybridStructure, low_part_size>(setHS, firstL, secondL, nTries, leftResult, r);
            // r.insert(r.end(), intersection_AND_HS.begin(), intersection_AND_HS.end());

            // vector<uint64_t> intersection_AND_HS_R; 
            // barbayKenyonHS<HybridStructure, low_part_size>(setHS, firstR, secondR, nTries, rightResult, r);
            // r.insert(r.end(), intersection_AND_HS_R.begin(), intersection_AND_HS_R.end());

            return true;       
        }
        existLeft = AND_High_BTrie_Log<HybridStructure, low_part_size>(setHS, nTries, maxLevel, nextLevel, leftNodes, leftResult, r);
        existRight = AND_High_BTrie_Log<HybridStructure, low_part_size>(setHS, nTries, maxLevel, nextLevel, rightNodes, rightResult, r);
    }
    if (existLeft || existRight)
        return true;
    return false;
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template <class HybridStructure, uint8_t low_part_size>
void Intersect_HS_High_BTrie_Log(vector<HybridStructure> &setHS, vector<uint64_t> &result, bool searchBF, bool runs, bool parallel){
    
    uint64_t max_level = setHS[0].getBinTrie()->getHeight();
    uint64_t n_tries = setHS.size();
    // max 16 relations
    bool activeHS[16] = { true, true, true, true,
                             true, true, true, true,
                             true, true, true, true,
                             true, true, true, true };
    uint64_t roots[16] = { 0, 0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0 };
    
    uint64_t curr_level = 0;
    
    // Hint of threads avaible in our CPU
    unsigned nb_threads_hint = std::thread::hardware_concurrency();
    uint64_t level_of_cut     = floor(log2(nb_threads_hint));
    unsigned int nb_threads   = pow(2, level_of_cut);
    
    vector<bool*>       activeTries2;
    vector<uint64_t*>   roots2;
    vector<uint64_t>    partial_solutions;
    
    result.reserve(1000000);
    
    if (parallel && level_of_cut <= max_level-2){
        uint64_t partial_int = 0;

        partialAND_HS(setHS, n_tries, max_level, curr_level, level_of_cut, roots, activeHS, partial_int, result,
                    partial_solutions, roots2, activeTries2);
        
        uint16_t real_threads = roots2.size();
        uint16_t init_threads = real_threads;
        // Go down one more level to fill all threads
        vector<bool*> next_active_tries;
        vector<uint64_t*> next_roots;
        vector<uint64_t> next_partial_solutions;
        uint64_t i = 0;
        
        while ((nb_threads - real_threads > 1) && (i < init_threads)) {
            partialAND_HS(setHS, n_tries, max_level, level_of_cut, level_of_cut+1, roots2[i], 
                        activeTries2[i], partial_solutions[i], result, 
                        next_partial_solutions, next_roots, next_active_tries);
            
            real_threads = init_threads + next_roots.size() - (++i);
        }

        if (real_threads > 0) {
            vector<uint64_t> next_init_level(next_roots.size(), level_of_cut+1);

            for(uint16_t j = i; j < init_threads; ++j) {
                next_roots.push_back(roots2[j]);
                // if (runs)
                next_active_tries.push_back(activeTries2[j]);
                next_partial_solutions.push_back(partial_solutions[j]);
                next_init_level.push_back(level_of_cut);
            }


        // Init the vector to contain threads solutions
            vector<vector<uint64_t>> threads_results(real_threads);
            for (auto tr: threads_results)
                tr.reserve(1000000);

            
            parallel_for(real_threads, real_threads, [&](int start, int end) {
                for (uint16_t threadId = start; threadId < end; ++threadId) {
                    
                        AND_High_BTrie_Log<HybridStructure, low_part_size>(setHS, n_tries, max_level, next_init_level[threadId], 
                            next_roots[threadId], next_partial_solutions[threadId], 
                            threads_results[threadId]);
                }
            });

            uint64_t output_size = result.size(); 
            vector<uint64_t> shifts(real_threads);
            uint64_t shift = result.size();
            for(uint64_t t = 0; t < real_threads; ++t){
                output_size += threads_results[t].size();
                shifts[t] = shift;
                shift += threads_results[t].size();
            }

            // Write in parallel threads result
            if (output_size > 450000){
                result.resize(output_size);
                parallel_for(real_threads, real_threads, [&](int start, int end) {
                    for (uint16_t threadId = start; threadId < end; ++threadId) {
                        for (uint64_t j = 0; j < threads_results[threadId].size(); ++j) {
                            result[j+shifts[threadId]] = threads_results[threadId][j];
                        } 
                    }        
                });
            } 
            else {
                // Concatenate solutions of threads
                for(uint64_t t=0; t < real_threads; ++t){
                    result.insert(  result.end(), 
                                    threads_results[t].begin(),
                                    threads_results[t].end()
                                );
                }

            }
        }
        // Free memory
        for (uint64_t i = 0; i < real_threads; ++i) {
            delete[] next_active_tries[i];
            delete[] next_roots[i];
        }
    }   else {        
        AND_High_BTrie_Log<HybridStructure, low_part_size>(setHS, n_tries, max_level, 0, roots, 0, result);
    }
    return;
}
