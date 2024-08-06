#ifndef INTERSECTION_BRUTE_FORCE_CPP
#define INTERSECTION_BRUTE_FORCE_CPP

#include <iostream> 
#include <vector>
#include <sdsl/int_vector.hpp>
#include "structures.hpp"


using namespace std;
using namespace sdsl;

template <uint8_t low_part_size>
uint64_t bruteForceSearchHS(sdsl::int_vector<low_part_size> *lowPart, uint64_t n, uint64_t x, uint64_t initial_position) {    
    if ((*lowPart)[initial_position] >= x)
        return initial_position;
    uint64_t value = (*lowPart)[n];
    if (value <= x)
        return n;
        
    uint64_t i = initial_position + 1; 
    value = (*lowPart)[i];
    while (i < n  && value <= x) {
        if(value == x)
            return i; 
        else
            i++;
    value = (*lowPart)[i];
    }

    return i; 
}

template <class HybridStructure, uint8_t low_part_size>
void IntersectBruteForce_HS(vector<HybridStructure> &setHS, uint64_t initial_position[], uint64_t final_position[], uint16_t k,
                            uint64_t high_part, vector<uint64_t> &intersection){
    uint64_t e = (*setHS[0].getLowPart())[initial_position[0]];
    uint8_t i = 1;
    uint8_t occr = 1;
    uint64_t size = final_position[i];

    while ( e != -1 ){
        sdsl::int_vector<low_part_size> *lowPart = setHS[i].getLowPart();
        uint64_t pos = bruteForceSearchHS<low_part_size>( lowPart, size, e, initial_position[i]);
        uint64_t value = (*lowPart)[pos];
        if (value == e ) {
            occr ++;
            initial_position[i] = pos;
            
            if(occr == k){
                intersection.push_back((high_part << low_part_size) | e);
            }
        } 
        
        if(occr == k || value != e ){
            if (pos == size && value <= e){
                e = -1;
                return;
            }
            if (occr == k){
                e = (*lowPart)[pos+1];
                initial_position[i] = pos + 1;
            }
            else{
                e = value;
                initial_position[i] = pos;
                
            }
            occr = 1;
        }
        i = (i+1)%k;
        size = final_position[i];
    }
    return;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------

template <class HybridStructure, uint8_t low_part_size>
void IntersectBruteForceRuns_HS(vector<HybridStructure> &setHS, uint64_t initial_position[], uint64_t final_position[], uint16_t k,
                            uint64_t high_part, vector<uint64_t> &intersection, uint64_t activePosHS[]){
    uint64_t  i = 0;
    uint64_t e = (*setHS[activePosHS[i]].getLowPart())[initial_position[activePosHS[i]]];
    uint64_t occr = 1;
    i++;
    uint64_t size = final_position[activePosHS[i]];

    while ( e != -1){
        uint64_t  posHS = activePosHS[i];
        sdsl::int_vector<low_part_size> *lowPart = setHS[posHS].getLowPart();
        uint64_t pos = bruteForceSearchHS<low_part_size>( lowPart, size, e, initial_position[posHS]);
        uint64_t value = (*lowPart)[pos];
        if (value == e ) {
            occr ++;
            initial_position[posHS] = pos;
            
            if(occr == k){
                intersection.push_back((high_part << low_part_size) | e);
            }
        }        
        if(occr == k || value != e ){
            if (pos == size && value <= e){
                e = -1;
                return;
            }
            if (occr == k){
                e = (*lowPart)[pos+1];
                initial_position[posHS] = pos + 1;
            }
            else{
                e = value;
                initial_position[posHS] = pos;
                
            }
            occr = 1;
        }
        i = (i+1)%k;
        size = final_position[activePosHS[i]];
    }
    return;
}

#endif