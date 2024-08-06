#ifndef BINARYSEARCH_HPP
#define BINARYSEARCH_HPP

#include <iostream>
#include <vector>
#include <algorithm>

inline int binarySearch(std::vector<uint64_t> * vec, uint64_t target) {
    int left = 0;
    int right = vec->size() - 1;

    while (left <= right) {
        int mid = left + (right - left) / 2;
        if ((*vec)[mid] == target) {
            return mid; 
        } else if ((*vec)[mid] < target) {            
            left = mid + 1;
        } else {
            right = mid - 1;
        }

    }
    return -1;
}

#endif