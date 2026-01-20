#include <cstdint>
#include <vector>
#include <sdsl/int_vector.hpp>

#ifndef STRUCTURES_H
#define STRUCTURES_H


using namespace std;
using namespace sdsl;

template <uint8_t low_part_size>
inline void printIntVector(sdsl::int_vector<low_part_size> *lowPart, int init, int end) {
    for(int i = init; i<end; i++){
        uint64_t value = (*lowPart)[i];
        cout << value;
        if(i < end-1)
            cout << " , ";
        else 
            cout << endl;
    }
}
inline void printVector(vector<uint64_t> &arr) {
    for(int i = 0; i < arr.size(); i++){
        cout << arr[i];
        if(i < arr.size()-1)
            cout << " -- ";
        else 
            cout << endl;
    }
}



// Barbay and Kenyon search data
struct search_data {
    uint64_t time;
    uint64_t total_time;
    double general_total_time;
    uint64_t longer_time;
    uint64_t total_longer_time; 
    uint64_t start_longer_time; 
    uint64_t end_longer_time; 
    uint64_t shorter_time;
    uint64_t total_shorter_time;  
    uint64_t start_shorter_time; 
    uint64_t end_shorter_time; 
    uint64_t number_searches;
};


#endif