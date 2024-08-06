#include <cstdint>

#ifndef STRUCTURES_H
#define STRUCTURES_H

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