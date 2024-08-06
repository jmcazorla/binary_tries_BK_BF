#ifndef BARBAY_AND_KENYON_H_
#define BARBAY_AND_KENYON_H_

#include <iostream> 
#include <vector>
#include <sdsl/int_vector.hpp>
#include "structures.hpp"

using namespace std;
using namespace sdsl;

void deltaBarbayKenyon(vector<vector<uint64_t>> &sets, uint16_t k, uint64_t &delta){
    vector<uint64_t> positions(k, 0);
    vector<uint64_t> final_positions(k, 0);
    vector<bool> range_covered_by(k, 0);
    int64_t e, l_value, pos_k;
    uint64_t final_quantity = 0;
    uint64_t previous = 0;
    uint64_t change_quantity = 0;

    for (size_t j = 0; j < k; j++) {
        if(sets[j][0] < e) {
            e = sets[j][0];
            l_value = sets[j][1];
        }
        final_positions[j] = sets[j].size();
    }

    while ( e != -1 ) {
        pos_k = 0;
        while ( pos_k < k ) {
            if (positions[pos_k] < final_positions[pos_k]) {                
                if(sets[pos_k][positions[pos_k]] == e) {
                    if(!range_covered_by[pos_k]) {
                        range_covered_by[pos_k] = 1;                        
                        change_quantity++;
                    }
                    positions[pos_k]++;
                    if (positions[pos_k] < final_positions[pos_k]) {
                        if ( sets[pos_k][positions[pos_k]] < l_value || l_value == -1) {
                            l_value = sets[pos_k][positions[pos_k]];
                        }                        
                    } else {
                        final_quantity++;
                    }
                } else if ( sets[pos_k][positions[pos_k]] < l_value || l_value == -1) {
                    l_value = sets[pos_k][positions[pos_k]];
                }
            } 
            pos_k++;
        }
        if(change_quantity == k) {
            delta += previous + 1;
            change_quantity=0;
            range_covered_by.assign(k, 0);
            previous = 0;
        } else if (previous == 0) {
            previous++;
        }
        if (final_quantity == k) {
            e = -1;
        } else {
            e = l_value;
            l_value = -1;
        }
    }
    if (previous == 1) delta++;
    return;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------

uint64_t binarySearch(vector<uint64_t> &array, uint64_t low, uint64_t high, uint64_t x, uint64_t size) {
	while (low <= high)
	{
        int mid = (high + low)/2;
		if (array[mid] == x)
			return mid;
		if (array[mid] < x)
			low = mid + 1;
		else
			high = mid - 1;
	}
    if (low > size - 1){
        return high;
    }
    if (high < 0){
        return low;
    }

    return low;
}

uint64_t exponentialSearch(vector<uint64_t> arr, uint64_t n, uint64_t x, uint64_t initial_position) { 
    if (arr[initial_position] >= x) 
        return initial_position; 
    if (arr[n-1] <= x)
        return n-1;
  
    uint64_t i = initial_position + 1; 
    while (i < n && arr[i] <= x) 
        i = i*2; 
    return binarySearch(arr, i/2, min(i, n-1), x, n); 
}

void barbayKenyon(vector<vector<uint64_t>> &sets, uint16_t k, vector<uint64_t> &intersection){
    vector<uint64_t> positions(k, 0);
    // eliminator element in [0,0], first element of first set
    uint64_t e = sets[0][0];
    // Set index
    uint64_t i = 1;
    // ocurrences of e
    uint64_t occr = 1;
    // Init actual elements and size of initial set
    uint64_t size = sets[i].size();
    while ( e != -1 ){
        // position of e in ith-set
        uint64_t pos = exponentialSearch(sets[i], size, e, positions[i]);
        if (sets[i][pos] == e){
            occr ++;
            positions[i] = pos;            
            if(occr == k){
                intersection.push_back(e);
            }
        } 
        if(occr == k || sets[i][pos] != e){
            // Position remain and size of next set
            uint64_t next_set_pos = positions[i];
            uint64_t next_set_size = sets[i].size();
            // No elements remain in the smallest set
            if (pos == next_set_size-1){          
                e = -1;
                return;
            }
            // e it's part of intersection      
            if (occr == k){
                e = sets[i][pos+1];
                positions[i] = pos + 1;
            }
            // e is not found in actual set
            else{
                // pos it's a succesor index of e 
                e = sets[i][pos];
                positions[i] = pos;                
            }
            // restart occurrences
            occr = 1;
        }
        // Cyclical index of sets
        i = (i+1)%k;
        size = sets[i].size();
    }
    return;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------

template <uint8_t low_part_size>
uint64_t binarySearchHS(sdsl::int_vector<low_part_size> *lowPart, uint64_t low, uint64_t high, uint64_t x, uint64_t size) {
    while (low <= high)
	{
        int mid = (high + low)/2;
        uint64_t value = (*lowPart)[mid];
		if (value == x)
			return mid;
		if (value < x)
			low = mid + 1;
		else
			high = mid - 1;
	}
    if (low > size){
        return high;
    }
    if (high < 0){
        return low;
    }    
    return low;
}

template <uint8_t low_part_size>
uint64_t exponentialSearchHS(sdsl::int_vector<low_part_size> *lowPart, uint64_t n, uint64_t x, uint64_t initial_position) { 
    if ((*lowPart)[initial_position] >= x) 
        return initial_position;
    if ((*lowPart)[n] <= x)
        return n;
    uint64_t i = 1; 
    uint64_t end = n - initial_position;
    while (i <= end && (*lowPart)[initial_position + i] <= x)
        i = i*2; 
    uint64_t low = initial_position + i/2;
    return binarySearchHS(lowPart, low, min(initial_position + i, n), x, n + 1); 
}

template <class HybridStructure, uint8_t low_part_size>
void barbayKenyonHS(vector<HybridStructure> &setHS, uint64_t initial_position[], uint64_t final_position[], uint16_t k,
                    uint64_t high_part, vector<uint64_t> &intersection){
    uint64_t e = (*setHS[0].getLowPart())[initial_position[0]];
    // Set index
    uint64_t i = 1;
    // ocurrences of e
    uint64_t occr = 1;
    // Init actual elements and size of initial set
    uint64_t size = final_position[i];
    while ( e != -1 ){    
        sdsl::int_vector<low_part_size> *lowPart = setHS[i].getLowPart();
        // position of e in ith-set
        uint64_t pos = exponentialSearchHS<low_part_size>( lowPart, size, e, initial_position[i]);
        uint64_t value = (*lowPart)[pos];
        if (value == e ) {
            occr ++;
            initial_position[i] = pos;            
            if(occr == k){
                intersection.push_back((high_part << low_part_size) | e);
            }
        }        
        if(occr == k || value != e ){
            if (pos >= size){
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

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------

template <class HybridStructure, uint8_t low_part_size>
void barbayKenyonRunsHS(vector<HybridStructure> &setHS, uint64_t initial_position[], uint64_t final_position[], uint16_t k,
                    uint64_t high_part, vector<uint64_t> &intersection, uint64_t activePosHS[]){
    uint64_t  i = 0;
    uint64_t e = (*setHS[activePosHS[i]].getLowPart())[initial_position[activePosHS[i]]];
    uint64_t occr = 1;
    i++;
    uint64_t size = final_position[activePosHS[i]];

    while ( e != -1 ){
        uint64_t  posHS = activePosHS[i];
        sdsl::int_vector<low_part_size> *lowPart = setHS[posHS].getLowPart();
        uint64_t pos = exponentialSearchHS<low_part_size>( lowPart, size, e, initial_position[posHS]);
        uint64_t value = (*lowPart)[pos];
        if (value == e ) {
            occr ++;
            initial_position[posHS] = pos;            
            if(occr == k){
                intersection.push_back((high_part << low_part_size) | e);
            }
        }        
        if(occr == k || value != e ){
            if (pos >= size){
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