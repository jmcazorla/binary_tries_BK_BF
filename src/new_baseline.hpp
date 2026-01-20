#include <iostream>
#include <vector>
#include <algorithm>
#include <vector>
#include <cstdint>

/**
 * binary search.
 */
size_t binary_search(const std::vector<uint64_t>& vec, size_t low, uint64_t value) {
    size_t high = vec.size();
    while (low < high) {
        size_t mid = low + (high - low) / 2;
        if (vec[mid] < value) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low;
}

/**
 * Binary search within a specific range.
 */
size_t binary_search_rank(const std::vector<uint64_t>& vec, size_t low, size_t high, uint64_t value) {
    while (low < high) {
        size_t mid = low + (high - low) / 2;
        if (vec[mid] < value) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }
    return low;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * Swapping SvS intersection of two sets.
 */
void intersect_swapping(const std::vector<uint64_t>& set_a, const std::vector<uint64_t>& set_b,
                                        std::vector<uint64_t>& result) {
    size_t idx_a = 0;
    size_t idx_b = 0;
    size_t size_a = set_a.size();
    size_t size_b = set_b.size();

    while (idx_a < size_a && idx_b < size_b) {
        if (set_a[idx_a] == set_b[idx_b]) {
            result.push_back(set_a[idx_a]);
            idx_a++;
            idx_b++;
        } 
        else {
            size_t remaining_a = size_a - idx_a;
            size_t remaining_b = size_b - idx_b;
            
            // Basic pointer advancement
            if(set_a[idx_a] > set_b[idx_b])
                ++idx_b;
            else
                ++idx_a;
            // Swapping logic: Pick the "eliminator" from the set with fewer remaining elements
            if (remaining_a <= remaining_b) {
                idx_b = binary_search(set_b, idx_b, set_a[idx_a]);
            } else  {
                idx_a = binary_search(set_a, idx_a, set_b[idx_b]);
            } 
        }
    }
}

/**
 * Swapping SvS algorithm for k sets.
 * Intersects sets two-by-two starting from the smallest ones.
 */
std::vector<uint64_t> swapping_svs(std::vector<std::vector<uint64_t>>& sets) {
    if (sets.empty()) return {};
    // Sort sets by initial size
    std::sort(sets.begin(), sets.end(), [](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
        return a.size() < b.size();
    });

    std::vector<uint64_t> candidate = sets[0];

    for (size_t i = 1; i < sets.size(); ++i) {
        if (candidate.empty()) break;

        std::vector<uint64_t> next_candidate;
        intersect_swapping(candidate, sets[i], next_candidate);
        candidate = std::move(next_candidate);
    }

    return candidate;
}

//----------------------------------------------------------------------------------------------------------------------------------------------------

/**
 * Recursive So_BaezaYates intersection between two sets.
 */
void BY_intersect_sorted(
    const std::vector<uint64_t>& setA, size_t minA, size_t maxA,
    const std::vector<uint64_t>& setB, size_t minB, size_t maxB,
    std::vector<uint64_t>& result) 
{
    if (minA > maxA || minB > maxB) {
        return;
    }

    size_t m = minA + (maxA - minA) / 2;
    uint64_t medianA = setA[m];
    // Find the rank of medianA in set B
    size_t r = binary_search_rank(setB, minB, maxB + 1, medianA);
    //left sub-problem
    if (m > 0 && r > 0) {
        BY_intersect_sorted(setA, minA, m - 1, setB, minB, r - 1, result);
    }

    if (r <= maxB && setB[r] == medianA) {
        result.push_back(medianA);
    }
    //right sub-problem
    BY_intersect_sorted(setA, m + 1, maxA, setB, r, maxB, result);
}

/**
 * So_BaezaYates algorithm for k sets.
 * Variant of Baeza-Yates that ensures sorted output without explicit sorting of intermediate results.
 */
std::vector<uint64_t> so_baeza_yates(std::vector<std::vector<uint64_t>>& sets) {
    if (sets.empty()) return {};

    std::sort(sets.begin(), sets.end(), [](const std::vector<uint64_t>& a, const std::vector<uint64_t>& b) {
        return a.size() < b.size();
    });

    std::vector<uint64_t> candidate = sets[0];

    for (size_t i = 1; i < sets.size(); ++i) {
        if (candidate.empty()) break;

        std::vector<uint64_t> next_candidate;
        BY_intersect_sorted(
            candidate, 0, candidate.size() - 1,
            sets[i], 0, sets[i].size() - 1,
            next_candidate
        );    
        candidate = std::move(next_candidate);
    }
    return candidate;
}



//----------------------------------------------------------------------------------------------------------------------------------------------------



/**
 * Galloping Search (Exponential Search).
 */
size_t galloping_search(const std::vector<uint64_t>& vec, size_t low, uint64_t value) {
    size_t n = vec.size();
    if (low >= n) return n;
    if (vec[low] >= value) return low;

    size_t step = 1;
    size_t high = low + step;

    // Expansion phase (galoping)
    while (high < n && vec[high] < value) {
        low = high;
        step <<= 1;
        high = low + step;
    }

    if (high >= n) high = n - 1;

    // Binary search phase in the bounded range [low, high]
    while (low <= high) {
        size_t mid = low + (high - low) / 2;
        if (vec[mid] < value) {
            low = mid + 1;
        } else {
            if (mid == 0 || vec[mid - 1] < value) return mid;
            high = mid - 1;
        }
    }
    return n;
}

/**
 * Small_Adaptive algorithm for k sorted sets.
 */
std::vector<uint64_t> small_adaptive(std::vector<std::vector<uint64_t>>& sets) {
    size_t k = sets.size();
    if (k == 0) return {};

    std::vector<uint64_t> result;
    std::vector<size_t> current(k, 0);
    std::vector<size_t> order(k);
    for (size_t i = 0; i < k; ++i) order[i] = i;

    while (true) {
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return (sets[a].size() - current[a]) < (sets[b].size() - current[b]);
        });

        size_t s_idx = order[0];
        if (current[s_idx] >= sets[s_idx].size()) break;

        uint64_t e = sets[s_idx][current[s_idx]];
        current[s_idx]++; 

        // Search for 'e' in the remaining k-1 sets
        size_t idx = 1;
        for (; idx < k; ++idx) {
            size_t curr_set_idx = order[idx];
            
            current[curr_set_idx] = galloping_search(sets[curr_set_idx], current[curr_set_idx], e);
            
            if (current[curr_set_idx] >= sets[curr_set_idx].size() || sets[curr_set_idx][current[curr_set_idx]] != e) {
                break; 
            }
        }

        if (idx == k) {
            result.push_back(e);
            for (size_t i = 1; i < k; ++i) {
                current[order[i]]++;
            }
        }
    }

    return result;
}