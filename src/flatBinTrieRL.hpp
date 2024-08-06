#ifndef FLAT_BINTRIE_BV_RL
#define FLAT_BINTRIE_BV_RL

#include <iostream>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <vector>
#include <queue>
#include <math.h>
#include <thread>
#include "binarySearch.hpp"

using namespace sdsl;
using namespace std;

template <class rankType = rank_support_v<1>>
class flatBinTrieRL{
    private:
        uint16_t height; // original height of trie
        uint16_t height_with_runs; // height with runs encoded
        bool empty_trie = false;
        bool runs_encoded;

        sdsl::bit_vector *bTrie;
        sdsl::bit_vector *lastLevel;     
        rankType b_rank;
        uint64_t* level_pos;
        rankType lL_rank;
        

    public:

        flatBinTrieRL() = default;

        ~flatBinTrieRL() = default;

        flatBinTrieRL(vector<uint64_t> &set, uint64_t u) {
            flatBinTrieRL::runs_encoded = false;
            uint32_t n = set.size();

            uint16_t height = floor(log2(u - 1)) +  1;
            flatBinTrieRL::height = height;
            flatBinTrieRL::level_pos = new uint64_t[height];
            
            uint64_t max_nodes            = 2 * (pow(2, height) - 1);
            uint64_t max_nodes_last_level = pow(2, height);  
            
            flatBinTrieRL::bTrie     = new bit_vector(max_nodes, 0);
            flatBinTrieRL::lastLevel = new bit_vector(max_nodes_last_level, 0);
            
            queue<tuple<uint64_t, uint64_t, uint64_t>> q;         
            // add all set to split
            tuple<uint64_t, uint64_t, uint64_t> split(0, n-1, n);
            q.push(split);

            uint16_t level            = 0;
            uint64_t nodes_curr_level = 1; 
            uint64_t count_nodes      = 0;
            uint64_t nodes_next_level = 0;
            uint64_t index            = 0;
            uint64_t total_nodes      = 0;
            uint64_t nodes_last_level = 0;

            while (!q.empty()) {
                count_nodes++; // count node visited
                // Get node to write
                tuple<uint64_t, uint64_t, uint64_t> s = q.front();
                q.pop(); 
                
                uint64_t l, r, n;
                std::tie(l, r, n) = s;
                uint64_t left_elements  = 0;
                uint64_t right_elements = 0;

                // j-th most significative bit
                uint8_t j = height - level;
                uint64_t ll, lr, rl, rr;
                for (uint64_t i = l; i < r+1; ++i) {
                    if ((set[i] >> j-1) & 1) {                        
                        right_elements = r-i+1;
                        rl = i;
                        rr = r;
                        break;
                    }
                    else {
                        if (i == l){
                            ll = l;
                        }
                        lr = i;    
                        left_elements++;
                    }
                }
                // Add to queue split sets and write nodes
                tuple<uint64_t,uint64_t,uint64_t> left_split;
                tuple<uint64_t,uint64_t,uint64_t> right_split;
                // left child
                if (left_elements > 0) {
                    // write 1
                    if (level == height-1)
                        (*lastLevel)[index] = 1;
                    else
                        (*bTrie)[index] = 1;
                    tuple<uint64_t,uint64_t,uint64_t> left_split(ll, lr, left_elements);
                    q.push(left_split);
                    nodes_next_level++;
                    index++;
                    total_nodes++;
                }
                else {
                    index++;
                    total_nodes++;
                }
                // right child
                if (right_elements > 0) {
                    // write 1
                    if (level == height-1)
                        (*lastLevel)[index] = 1;
                    else
                        (*bTrie)[index] = 1;
                    tuple<uint64_t,uint64_t,uint64_t> right_split(rl, rr, right_elements);
                    q.push(right_split);
                    nodes_next_level++;
                    index++;
                    total_nodes++;
                }
                else {
                    index++;
                    total_nodes++;
                }

                if (count_nodes == nodes_curr_level) {
                    // bTrie[level].resize(2*count_nodes);'
                    level_pos[level] = index;
                    if (level == height-2){
                        nodes_last_level = nodes_next_level;
                        index = 0;
                    }
                    
                    nodes_curr_level = nodes_next_level;
                    nodes_next_level = 0;
                    count_nodes = 0;
                    level++;
                    
                }
                if (level == flatBinTrieRL::height) {
                    break;
                }
            }

            flatBinTrieRL::bTrie -> resize(total_nodes - 2*nodes_last_level);
            flatBinTrieRL::lastLevel -> resize(2*nodes_last_level);
            flatBinTrieRL::b_rank = rankType(bTrie);
            flatBinTrieRL::lL_rank = rankType(lastLevel);
        }

        
        flatBinTrieRL(int_vector<> &set) {
            flatBinTrieRL::runs_encoded = false;
            uint32_t n = set.size();
            
            util::bit_compress(set);
            flatBinTrieRL::height = (uint16_t)set.width();
            flatBinTrieRL::level_pos = new uint64_t[height];
            
            uint64_t max_nodes = 2*(pow(2, height+1) - 1);
            flatBinTrieRL::bTrie = new bit_vector(max_nodes, 0); 
            

            queue<tuple<uint64_t, uint64_t, uint64_t>> q;
            
            // add all set to split
            tuple<uint64_t, uint64_t, uint64_t> split(0, n-1, n);
            q.push(split);

            uint16_t level            = 0;
            uint64_t nodes_curr_level = 1; 
            uint64_t count_nodes      = 0;
            uint64_t nodes_next_level = 0;
            uint64_t index            = 0;
            uint64_t total_nodes      = 0;

            while (!q.empty()) {
                count_nodes++; // count node visited
                // Get node to write
                tuple<uint64_t, uint64_t, uint64_t> s = q.front();
                q.pop(); 
                
                uint64_t l, r, n;
                std::tie(l, r, n) = s;
                uint64_t left_elements  = 0;
                uint64_t right_elements = 0;

                // j-th most significative bit
                uint8_t j = height - level;
                uint64_t ll, lr, rl, rr;
                for (uint64_t i = l; i < r+1; ++i) {
                    if ((set[i] >> j-1) & 1) {                        
                        right_elements = r-i+1;
                        rl = i;
                        rr = r;
                        break;
                    }
                    else {
                        if (i == l){
                            ll = l;
                        }
                        lr = i;    
                        left_elements++;
                    }
                }
                // Add to queue split sets and write nodes
                tuple<uint64_t,uint64_t,uint64_t> left_split;
                tuple<uint64_t,uint64_t,uint64_t> right_split;
                // left child
                if (left_elements > 0) {
                    // write 1
                    (*bTrie)[index] = 1;
                    tuple<uint64_t,uint64_t,uint64_t> left_split(ll, lr, left_elements);
                    q.push(left_split);
                    nodes_next_level++;
                    index++;
                    total_nodes++;
                }
                else {
                    // write 0
                    (*bTrie)[index] = 0;
                    index++;
                    total_nodes++;
                }
                // right child
                if (right_elements > 0) {
                    // write 1
                    (*bTrie)[index] = 1;
                    tuple<uint64_t,uint64_t,uint64_t> right_split(rl, rr, right_elements);
                    q.push(right_split);
                    nodes_next_level++;
                    index++;
                    total_nodes++;
                }
                else {
                    // write 0
                    (*bTrie)[index] = 0;
                    index++;
                    total_nodes++;
                }

                if (count_nodes == nodes_curr_level) {
                    level_pos[level] = index;
                    nodes_curr_level = nodes_next_level;
                    nodes_next_level = 0;
                    count_nodes = 0;
                    level++;
                    
                }
                if (level == flatBinTrieRL::height) {
                    break;
                }
            }
            flatBinTrieRL::bTrie -> resize(total_nodes);
            flatBinTrieRL::b_rank = rankType(bTrie);
            flatBinTrieRL::lL_rank = rankType(lastLevel);
        };


        flatBinTrieRL(vector<uint64_t> ones_to_write[], uint16_t height, 
                    uint64_t* level_pos, 
                    bool runs_encoded) {
            flatBinTrieRL::runs_encoded = runs_encoded;
            flatBinTrieRL::height = height;
            
            uint64_t bits_n = 0;
            uint16_t all_levels_empty = 0;
            uint16_t max_level_not_empty = 0;
            
            for (uint16_t level = 0; level < height; ++level) {
                bits_n += level_pos[level];
                if (level_pos[level] > 0){
                    max_level_not_empty = level;
                }
                else {
                    all_levels_empty++;
                }
                    
            }

            if (all_levels_empty == flatBinTrieRL::height) {
                flatBinTrieRL::empty_trie = true;
            }
            else {
                flatBinTrieRL::empty_trie = false;
            }

            flatBinTrieRL::height_with_runs = max_level_not_empty + 1;

            flatBinTrieRL::bTrie     = new bit_vector(bits_n - level_pos[max_level_not_empty], 0);
            flatBinTrieRL::lastLevel = new bit_vector(level_pos[max_level_not_empty], 0); 
            flatBinTrieRL::level_pos = new uint64_t[height];

            uint64_t global_level_pos = 0;
            for (uint16_t level = 0; level < height; ++level) {
                for (uint64_t i = 0; i < ones_to_write[level].size(); ++i) {
                    if (level == max_level_not_empty) {
                        uint64_t pos = ones_to_write[level][i];
                        (*flatBinTrieRL::lastLevel)[pos] = 1;
                    }
                    else {
                        uint64_t pos = global_level_pos + ones_to_write[level][i];
                        (*flatBinTrieRL::bTrie)[pos] = 1;
                    }
                }
                // cout << endl;
                if (level == max_level_not_empty) {
                    flatBinTrieRL::level_pos[level] = level_pos[level];
                }
                else {
                    global_level_pos += level_pos[level];
                    flatBinTrieRL::level_pos[level] = global_level_pos;
                }
                
            }
            flatBinTrieRL::b_rank = rankType(flatBinTrieRL::bTrie);
            flatBinTrieRL::lL_rank = rankType(flatBinTrieRL::lastLevel);
        };
        

        // free memory of bit vectors
        inline void free(){
            delete flatBinTrieRL::bTrie;
            delete flatBinTrieRL::lastLevel;
        }

        // return number of elements of bit_vector
        inline uint64_t size(){
            return bTrie->size();
        };


        inline uint16_t getHeight(){
            return flatBinTrieRL::height;
        };

        // Return number of elements coded in a trie
        inline uint64_t elements_coded() {
            uint64_t counter_ones = 0;
            for (uint64_t i = 0; i < flatBinTrieRL::lastLevel -> size(); ++i) {
                if ((*flatBinTrieRL::lastLevel)[i] == 1) 
                    counter_ones++;
            }
            return counter_ones;
        }

        
        inline uint64_t getNode(uint64_t &node_id) {
            uint64_t node = 0;
            uint64_t pos;
            if ((2*node_id) >= (flatBinTrieRL::bTrie -> size())) {
                pos = (2*node_id) - (flatBinTrieRL::bTrie -> size());
                if ((*lastLevel)[pos])
                    node = (node | (1ULL << 1));

                if ((*lastLevel)[pos + 1])
                    node = (node | (1ULL << 0));
            }
            else { 
                pos = 2*node_id;
                if ((*bTrie)[pos])
                    node = (node | (1ULL << 1));

                if ((*bTrie)[pos + 1])
                    node = (node | (1ULL << 0));
            }

            return node;
        };

        // GetNode in h-1 levels
        inline uint64_t getNode1(uint64_t &node_id) {
                return (((*bTrie)[2 * node_id]) << 1) | (*bTrie)[(2 * node_id)+1];
        }

        // GetNode in last level
        inline int64_t getNode2(uint64_t &node_id) {
            return ((*lastLevel)[2*node_id -flatBinTrieRL::bTrie -> size()] << 1) | (*lastLevel)[(2*node_id -flatBinTrieRL::bTrie -> size())+1];
        }


        inline uint64_t getLeftChild(uint64_t &node_id) {
                return flatBinTrieRL::b_rank((2*node_id) + 1);;
        };


        inline uint64_t getRightChild(uint64_t &node_id) {
                return flatBinTrieRL::b_rank((2*node_id) + 2);
        };

        inline uint64_t getChildInLastLevel(uint64_t &node_id) {
                uint64_t pos = 2*node_id - (bTrie -> size());
                return flatBinTrieRL::lL_rank(pos + 1);
        };

        inline uint64_t size_in_bytes() {
            uint64_t bv_size = sdsl::size_in_bytes(*(flatBinTrieRL::bTrie));
            uint64_t lastL_size = sdsl::size_in_bytes(*(flatBinTrieRL::lastLevel));
            uint64_t rank_size = sdsl::size_in_bytes(flatBinTrieRL::b_rank);
            uint64_t lL_rank_size  = sdsl::size_in_bytes(flatBinTrieRL::lL_rank);
            return bv_size +
                    rank_size +
                    lastL_size +
                    lL_rank_size +
                    2 * sizeof(bool) +
                    2 * sizeof(uint8_t);
        };

        // return size of bytes of all structure
        uint64_t serialize(std::ostream &out) {

            out.write(reinterpret_cast<char*>(&height)          , sizeof(height));
            out.write(reinterpret_cast<char*>(&height_with_runs), sizeof(height_with_runs));
            out.write(reinterpret_cast<char*>(&empty_trie)      , sizeof(empty_trie));
            out.write(reinterpret_cast<char*>(&runs_encoded)    , sizeof(runs_encoded));

            uint64_t bvs_size, rank_size;
            bvs_size  = bTrie -> serialize(out) + lastLevel -> serialize(out);
            rank_size = b_rank.serialize(out) + lL_rank.serialize(out);

            return bvs_size + rank_size + 
                   sizeof(height) + sizeof(height_with_runs) +
                   sizeof(empty_trie) + sizeof(runs_encoded);
        }

        // load structure from in stream
        void load(std::istream &in){
            in.read(reinterpret_cast<char*>(&height)          , sizeof(height));
            in.read(reinterpret_cast<char*>(&height_with_runs), sizeof(height_with_runs));
            in.read(reinterpret_cast<char*>(&empty_trie)      , sizeof(empty_trie));
            in.read(reinterpret_cast<char*>(&runs_encoded)    , sizeof(runs_encoded));

            flatBinTrieRL::bTrie = new sdsl::bit_vector();
            flatBinTrieRL::lastLevel = new sdsl::bit_vector();

            flatBinTrieRL::bTrie     -> load(in);
            flatBinTrieRL::lastLevel -> load(in);
            b_rank.load(in, bTrie);
            lL_rank.load(in, lastLevel);
        }


        inline void print() {
            uint64_t i = 0;
            for (uint16_t level=0; level < flatBinTrieRL::height; ++level) {
                uint64_t next_level_pos = flatBinTrieRL::level_pos[level];
                if (level == flatBinTrieRL::height - 1) {
                    i = 0;
                }
                while (i < next_level_pos) {
                    if (level < flatBinTrieRL::height - 1)
                        cout << (*flatBinTrieRL::bTrie)[i] << (*flatBinTrieRL::bTrie)[i+1] << " ";
                    else{ 
                        cout << (*flatBinTrieRL::lastLevel)[i] << (*flatBinTrieRL::lastLevel)[i+1] << " ";
                    }
                    ++(++i);
                }
                cout << endl;
            }
        };


        void writeCompressTrie(vector<uint64_t> ones_to_write[], vector<uint64_t> * ones_del_l_level,
                        uint64_t* level_pos,  std::vector<uint64_t> * posOnes,
                        uint16_t curr_level, uint64_t node_id, bool &its11){

            if (curr_level == (flatBinTrieRL::height-1)) {
                 uint64_t node = getNode2(node_id);
                uint64_t posL,posR;
                
                if (node == 0b11 ) {
                    posL = getChildInLastLevel(node_id);
                    posR = posL + 1;
                    int index = binarySearch(posOnes, posL); 
                    if (index != -1 && index < (posOnes->size() - 1))
                        its11 = (*posOnes)[index+1] == posR ? true : false ; 
                }
                if(its11){
                    ones_del_l_level->push_back(posL);
                    ones_del_l_level->push_back(posR);
                }
                if (node == 0b10 || (node == 0b11 && !its11)){
                    ones_to_write[curr_level].push_back(level_pos[curr_level]);
                }
                level_pos[curr_level] += 1;
                if (node == 0b01 || (node == 0b11 && !its11)) {
                    ones_to_write[curr_level].push_back(level_pos[curr_level]);
                }
                level_pos[curr_level] += 1;
                return;
            }
            uint64_t node = getNode1(node_id);
            uint16_t next_level = curr_level + 1;

            bool its11_l = false;
            bool its11_r = false;
            bool actualIts11 = false;

            if (node == 0b11) {
                actualIts11 = true;

                uint64_t l_child = getLeftChild(node_id);
                uint64_t r_child = l_child + 1;
                writeCompressTrie(ones_to_write, ones_del_l_level, level_pos, posOnes, next_level, l_child, its11_l);
                writeCompressTrie(ones_to_write, ones_del_l_level, level_pos, posOnes, next_level, r_child, its11_r);
                
                its11 = true && its11_l && its11_r;
                if (its11) {
                    level_pos[next_level] -= 4;
                }
                else {
                    ones_to_write[curr_level].push_back(level_pos[curr_level]);
                    ones_to_write[curr_level].push_back(level_pos[curr_level] + 1);      
                }                
                level_pos[curr_level] += 2;
            }

            else if (node == 0b10 || node == 0b01){
                if (node == 0b10){
                    uint64_t l_child = getLeftChild(node_id);
                    ones_to_write[curr_level].push_back(level_pos[curr_level]);
                    writeCompressTrie(ones_to_write, ones_del_l_level, level_pos, posOnes, next_level, l_child, its11_l);
                }
                level_pos[curr_level] += 1;
                if (node == 0b01) {
                    uint64_t r_child = getRightChild(node_id);
                    ones_to_write[curr_level].push_back(level_pos[curr_level]);
                    writeCompressTrie(ones_to_write, ones_del_l_level, level_pos, posOnes, next_level, r_child, its11_r);
                }
                level_pos[curr_level] += 1;
            }
        };

        
        // Method write ones in bit vector
               
        inline void writeOnes(vector<uint64_t> ones_to_write[], uint64_t* level_pos){
            flatBinTrieRL::runs_encoded = true;
            uint64_t bits_n = 0;
            uint16_t last_level = 0;
            uint64_t bits_before_last_level;
            for (uint16_t level = 0; level < flatBinTrieRL::height; ++level) {
                bits_n += level_pos[level];
                if (level_pos[level] > 0) {
                    last_level = level; 
                }               
            }

            flatBinTrieRL::height_with_runs = last_level + 1;
            delete flatBinTrieRL::bTrie;
            delete flatBinTrieRL::lastLevel;
            
            bit_vector* _bTrie     = new bit_vector(bits_n - level_pos[last_level], 0);
            bit_vector* _lastLevel = new bit_vector(level_pos[last_level], 0);
            flatBinTrieRL::level_pos = new uint64_t[height];

            uint64_t global_level_pos = 0;
            for (uint16_t level = 0; level < height_with_runs; ++level) {
                for (uint64_t i = 0; i < ones_to_write[level].size(); ++i) {
                    if (level == last_level){
                        uint64_t pos = ones_to_write[level][i];
                        (*_lastLevel)[pos] = 1;
                    }
                    else {
                        uint64_t pos = global_level_pos + ones_to_write[level][i];
                        (*_bTrie)[pos] = 1;
                    }
                }

                if (level == last_level) {
                    flatBinTrieRL::level_pos[level] = level_pos[level];
                }
                else {
                    global_level_pos += level_pos[level];
                    flatBinTrieRL::level_pos[level] = global_level_pos;
                }
            }


            flatBinTrieRL::bTrie = new bit_vector(*_bTrie);
            delete _bTrie; 
            flatBinTrieRL::lastLevel = new bit_vector(*_lastLevel);
            delete _lastLevel;
            flatBinTrieRL::b_rank = rankType(flatBinTrieRL::bTrie);
            flatBinTrieRL::lL_rank = rankType(flatBinTrieRL::lastLevel);
            
            delete[] level_pos;
        };


        inline vector<uint64_t> * encodeRuns(std::vector<uint64_t> * posOnes) {
            vector<uint64_t> ones_to_write[flatBinTrieRL::height];
            vector<uint64_t> * ones_del_l_level = new vector<uint64_t>();
            flatBinTrieRL::level_pos = new uint64_t[flatBinTrieRL::height];
            for(uint64_t i = 0; i < flatBinTrieRL::height; ++i) level_pos[i] = 0;

            bool itsOneOne = false;

            flatBinTrieRL::writeCompressTrie(ones_to_write, ones_del_l_level, level_pos, posOnes, 0, 0, itsOneOne);
            flatBinTrieRL::writeOnes(ones_to_write, level_pos);
            return ones_del_l_level;
        };

        
        inline void recursiveDecode(vector<uint64_t> &decoded, uint64_t partial_int, uint64_t node_id, uint16_t curr_level) {
            
            if (curr_level == flatBinTrieRL::height) {
                decoded.push_back(partial_int);
                return;
            }

            uint64_t node = getNode(node_id);
            uint16_t next_level = curr_level + 1;

            uint64_t leftResult = partial_int;
            uint64_t rightResult = partial_int;

            if (node == 0b10 || node == 0b11) {
                uint64_t left_child = getLeftChild(node_id);
                recursiveDecode(decoded, leftResult, left_child, next_level);
            }
            if (node == 0b01 || node == 0b11) {
                rightResult = (rightResult | (1ULL << (getHeight() - curr_level - 1)));
                uint64_t right_child = getRightChild(node_id);
                recursiveDecode(decoded, rightResult, right_child, next_level);
            }
        };


        inline void runsRecursiveDecode(vector<uint64_t> &decoded, uint64_t partial_int, uint64_t node_id, uint16_t curr_level) {
            if (curr_level == flatBinTrieRL::height) {
                decoded.push_back(partial_int);
                return;
            }

            uint64_t node = getNode(node_id);
            uint16_t next_level = curr_level + 1;

            if (node == 0b00) { 
                uint64_t below = partial_int;
                uint64_t range = ((uint64_t)1 << (getHeight() - curr_level)) - 1;
                uint64_t above = (partial_int | range);
                for (uint64_t i = below; i <= above; ++i) {
                    decoded.push_back(i);
                }
                return;
            }

            uint64_t leftResult = partial_int;
            uint64_t rightResult = partial_int;

            if (node == 0b10 || node == 0b11) {
                uint64_t left_child = getLeftChild(node_id);
                runsRecursiveDecode(decoded, leftResult, left_child, next_level);
            }
            if (node == 0b01 || node == 0b11) {
                rightResult = (rightResult | (1ULL << (getHeight() - curr_level - 1)));
                uint64_t right_child = getRightChild(node_id);
                runsRecursiveDecode(decoded, rightResult, right_child, next_level);
            }
        }


        inline void decode( vector<uint64_t> &decoded) {
            if (flatBinTrieRL::runs_encoded) {
                if (flatBinTrieRL::empty_trie) {
                    return;
                }
                else {
                    uint64_t partial_int = 0x00;
                    runsRecursiveDecode(decoded, partial_int, 0, 0);
                }
                
                
            }
            else {
                uint64_t partial_int = 0x00;
                recursiveDecode(decoded, partial_int, 0, 0);
            }
        }

        // If runs are encoded, this measure it's trie-run
        inline uint32_t trieMeasure() {
            uint64_t nEdgesLastLevel = 0;
            for (uint64_t i = 0; i < flatBinTrieRL::bTrie -> size(); ++i){
                if((*bTrie)[i] == 1) nEdgesLastLevel++;
            }
            for (uint64_t i = 0; i < flatBinTrieRL::lastLevel -> size(); ++i){
                if ((*lastLevel)[i] == 1) nEdgesLastLevel++;
            }
            return 
            nEdgesLastLevel;
        }
    
};

#endif
