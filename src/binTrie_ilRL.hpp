#ifndef BINTRIE_IL_RL
#define BINTRIE_IL_RL


#include <iostream>
#include <vector>
#include <queue>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/bit_vector_il.hpp>
#include "binarySearch.hpp"

using namespace std;
using namespace sdsl;


template <uint32_t block_size = 512>
class binTrie_ilRL{
    private:
        uint16_t height;
        uint16_t height_with_runs;
        bool empty_trie = false;
        bool runs_encoded;

        sdsl::bit_vector_il<block_size> *bTrie;
        sdsl::bit_vector_il<block_size> *lastLevel;     
        rank_support_il<1, block_size> b_rank;
        rank_support_il<1, block_size> lL_rank;
        uint64_t* level_pos;

    public:
        binTrie_ilRL() = default;

        ~binTrie_ilRL() = default;

        binTrie_ilRL(vector<uint64_t> &set, uint64_t u) {
            binTrie_ilRL::runs_encoded = false;
            uint32_t n = set.size();

            uint16_t height = floor(log2(u - 1)) +  1;
            binTrie_ilRL::height = height;
            binTrie_ilRL::level_pos = new uint64_t[height];
            
            uint64_t max_nodes            = 2 * (pow(2, height) - 1);
            uint64_t max_nodes_last_level = pow(2, height);  
            
            bit_vector* _bTrie     = new bit_vector(max_nodes, 0);
            bit_vector* _lastLevel = new bit_vector(max_nodes_last_level, 0);
            
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
                        (*_lastLevel)[index] = 1;
                    else
                        (*_bTrie)[index] = 1;
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
                        (*_lastLevel)[index] = 1;
                    else
                        (*_bTrie)[index] = 1;
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
                if (level == binTrie_ilRL::height) {
                    break;
                }
            }

            _bTrie -> resize(total_nodes - 2*nodes_last_level);
            _lastLevel -> resize(2*nodes_last_level);
            binTrie_ilRL::bTrie = new bit_vector_il<block_size>(*_bTrie);
            binTrie_ilRL::lastLevel = new bit_vector_il<block_size>(*_lastLevel);
            delete _bTrie;
            delete _lastLevel;
            binTrie_ilRL::b_rank = rank_support_il<1, block_size>(bTrie);
            binTrie_ilRL::lL_rank = rank_support_il<1, block_size>(lastLevel);
        }


        binTrie_ilRL(vector<uint64_t> ones_to_write[], uint16_t height, 
                    uint64_t* level_pos, 
                    bool runs_encoded) {
            binTrie_ilRL::runs_encoded = runs_encoded;
            binTrie_ilRL::height = height;
            
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

            if (all_levels_empty == binTrie_ilRL::height) {
                binTrie_ilRL::empty_trie = true;
            }
            else {
                binTrie_ilRL::empty_trie = false;
            }

            binTrie_ilRL::height_with_runs = max_level_not_empty + 1;

            bit_vector* _bTrie    = new bit_vector(bits_n - level_pos[max_level_not_empty], 0);
            bit_vector* _lastLevel = new bit_vector(level_pos[max_level_not_empty], 0); 
            binTrie_ilRL::level_pos = new uint64_t[height]; 

            uint64_t global_level_pos = 0;
            for (uint16_t level = 0; level < height; ++level) {
                for (uint64_t i = 0; i < ones_to_write[level].size(); ++i) {
                    if (level == max_level_not_empty) {
                        uint64_t pos = ones_to_write[level][i];
                        (*_lastLevel)[pos] = 1;
                    }
                    else {
                        uint64_t pos = global_level_pos + ones_to_write[level][i];
                        (*_bTrie)[pos] = 1;
                    }
                }
                if (level == max_level_not_empty) {
                    binTrie_ilRL::level_pos[level] = level_pos[level];
                }
                else {
                    global_level_pos += level_pos[level];
                    binTrie_ilRL::level_pos[level] = global_level_pos;
                }
                
            }
            binTrie_ilRL::bTrie = new bit_vector_il<block_size>(*_bTrie);
            binTrie_ilRL::lastLevel = new bit_vector_il<block_size>(*_lastLevel);
            delete _bTrie;
            delete _lastLevel;
            binTrie_ilRL::b_rank = rank_support_il<1, block_size>(binTrie_ilRL::bTrie);
            binTrie_ilRL::lL_rank = rank_support_il<1, block_size>(binTrie_ilRL::lastLevel);
        };

        // free memory of bit vectors
        inline void free(){
            delete binTrie_ilRL::bTrie;
            delete binTrie_ilRL::lastLevel;
        }

        // return number of elements of bit_vector
        inline uint64_t size(){
            return bTrie->size();
        };


        inline uint16_t getHeight(){
            return binTrie_ilRL::height;
        };

        // Return number of elements coded in a trie
        inline uint64_t elements_coded() {
            uint64_t counter_ones = 0;
            for (uint64_t i = 0; i < binTrie_ilRL::lastLevel -> size(); ++i) {
                if ((*binTrie_ilRL::lastLevel)[i] == 1) 
                    counter_ones++;
            }
            return counter_ones;
        }
        
        // Returns the node number up to the given level
        inline uint64_t elements_level(uint64_t level) {
            return binTrie_ilRL::level_pos[level];
        }

        inline uint64_t getNode(uint64_t node_id) {
            uint64_t node = 0;
            uint64_t pos;
            if ((2*node_id) >= (binTrie_ilRL::bTrie -> size())) {
                pos = (2*node_id) - (binTrie_ilRL::bTrie -> size());
                if ((*lastLevel)[pos])
                    node = (node | (1ULL << 1));

                if ((*lastLevel)[pos + 1] == 1)
                    node = (node | (1ULL << 0));
            }
            else { 
                pos = 2*node_id;
                if ((*bTrie)[pos])
                    node = (node | (1ULL << 1));

                if ((*bTrie)[pos + 1] == 1)
                    node = (node | (1ULL << 0));
            }

            return node;
        };

        inline uint64_t getNode1(uint64_t &node_id) {
            uint64_t pos  = 2 * node_id;
                return (((*bTrie)[pos]) << 1) | (*bTrie)[pos+1];
        }

        uint64_t getNode2(uint64_t &node_id) {
            uint64_t pos = 2*node_id - (bTrie -> size());
            return ((*lastLevel)[pos] << 1) | (*lastLevel)[pos+1];
        }


        inline uint64_t getLeftChild(uint64_t &node_id) {
                uint64_t rank = binTrie_ilRL::b_rank((2*node_id) + 1);
                return rank;
        };


        inline uint64_t getRightChild(uint64_t &node_id) {
                return binTrie_ilRL::b_rank((2*node_id) + 2);
        };

        inline uint64_t getChildInLastLevel(uint64_t &node_id) {
                uint64_t pos = 2*node_id - (bTrie -> size());
                return binTrie_ilRL::lL_rank(pos + 1);
        };

        inline uint64_t size_in_bytes() {
            uint64_t bv_size    = sdsl::size_in_bytes(*(binTrie_ilRL::bTrie));
            uint64_t lastL_size = sdsl::size_in_bytes(*(binTrie_ilRL::lastLevel));
            uint64_t rank_size  = sdsl::size_in_bytes(binTrie_ilRL::b_rank);
            uint64_t lL_rank_size  = sdsl::size_in_bytes(binTrie_ilRL::lL_rank);
            return bv_size +
                    rank_size +
                    lastL_size +
                    lL_rank_size +
                    2 * sizeof(bool) +
                    2 * sizeof(uint8_t);
        };


        // return size of bytes of all structure
        inline uint64_t serialize(std::ostream &out) {

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
        inline void load(std::istream &in){
            in.read(reinterpret_cast<char*>(&height)          , sizeof(height));
            in.read(reinterpret_cast<char*>(&height_with_runs), sizeof(height_with_runs));
            in.read(reinterpret_cast<char*>(&empty_trie)      , sizeof(empty_trie));
            in.read(reinterpret_cast<char*>(&runs_encoded)    , sizeof(runs_encoded));

            binTrie_ilRL::bTrie      = new sdsl::bit_vector_il<block_size>();
            binTrie_ilRL::lastLevel = new sdsl::bit_vector_il<block_size>();

            bTrie     -> load(in);
            lastLevel -> load(in);
            b_rank.load(in, bTrie);
            lL_rank.load(in, lastLevel);
        }


        inline void print() {
            uint64_t i = 0;
            for (uint16_t level=0; level < binTrie_ilRL::height; ++level) {
                uint64_t next_level_pos = binTrie_ilRL::level_pos[level];
                if (level == binTrie_ilRL::height - 1) {
                    i = 0;
                }
                while (i < next_level_pos) {
                    if (level < binTrie_ilRL::height - 1)
                        cout << (*binTrie_ilRL::bTrie)[i] << (*binTrie_ilRL::bTrie)[i+1] << " ";
                    else{ 
                        cout << (*binTrie_ilRL::lastLevel)[i] << (*binTrie_ilRL::lastLevel)[i+1] << " ";
                    }
                    ++(++i);
                }
                cout << endl;
            }
        };



        inline void writeCompressTrie(vector<uint64_t> ones_to_write[], vector<uint64_t> * ones_del_l_level,
                                uint64_t* aux_level_pos, std::vector<uint64_t> * posOnes, 
                                uint16_t curr_level, uint64_t node_id, bool &its11){   
            // End condition
            if (curr_level == (binTrie_ilRL::height-1)) {
                uint64_t node = getNode2(node_id);
                uint64_t posL,posR;
                
                if (node == 0b11 ) {
                    posL = getChildInLastLevel(node_id);
                    posR = posL + 1;
                    int index = binarySearch(posOnes, posL); 
                    if (index != -1 && index < (posOnes->size()-1))
                        its11 = (*posOnes)[index+1] == posR ? true : false ; 
                }
                if(its11){
                    ones_del_l_level->push_back(posL);
                    ones_del_l_level->push_back(posR);
                }
                if (node == 0b10 || (node == 0b11 && !its11)){
                    ones_to_write[curr_level].push_back(aux_level_pos[curr_level]);
                }
                aux_level_pos[curr_level] += 1;
                if (node == 0b01 || (node == 0b11 && !its11)) {
                    ones_to_write[curr_level].push_back(aux_level_pos[curr_level]);
                }
                aux_level_pos[curr_level] += 1;
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
                writeCompressTrie(ones_to_write, ones_del_l_level, aux_level_pos, posOnes, next_level, l_child, its11_l);
                writeCompressTrie(ones_to_write, ones_del_l_level, aux_level_pos, posOnes, next_level, r_child, its11_r);
                
                its11 = true && its11_l && its11_r;
                if (its11) {
                    aux_level_pos[next_level] -= 4;
                }
                else {
                    ones_to_write[curr_level].push_back(aux_level_pos[curr_level]);
                    ones_to_write[curr_level].push_back(aux_level_pos[curr_level] + 1);
                }
                aux_level_pos[curr_level] += 2;
            }

            else if (node == 0b10 || node == 0b01){
                if (node == 0b10){
                    uint64_t l_child = getLeftChild(node_id);
                    ones_to_write[curr_level].push_back(aux_level_pos[curr_level]);
                    writeCompressTrie(ones_to_write, ones_del_l_level, aux_level_pos, posOnes, next_level, l_child, its11_l);
                }
                aux_level_pos[curr_level] += 1;
                if (node == 0b01) {
                    uint64_t r_child = getRightChild(node_id);
                    ones_to_write[curr_level].push_back(aux_level_pos[curr_level]);
                    writeCompressTrie(ones_to_write, ones_del_l_level, aux_level_pos, posOnes, next_level, r_child, its11_r);
                }
                aux_level_pos[curr_level] += 1;
            }
        };

        // Method write ones in bit vector
        inline void writeOnes(vector<uint64_t> ones_to_write[], uint64_t* aux_level_pos){
            
            binTrie_ilRL::runs_encoded = true;
            uint64_t bits_n = 0;
            uint16_t last_level = 0;
            uint64_t bits_before_last_level;
            for (uint16_t level = 0; level < binTrie_ilRL::height; ++level) {
                bits_n += aux_level_pos[level];
                if (aux_level_pos[level] > 0) {
                    last_level = level; 
                }                
            }

            binTrie_ilRL::height_with_runs = last_level + 1;
            delete binTrie_ilRL::bTrie;
            delete binTrie_ilRL::lastLevel;
            
            bit_vector* _bTrie     = new bit_vector(bits_n - aux_level_pos[last_level], 0);
            bit_vector* _lastLevel = new bit_vector(level_pos[last_level], 0);
            binTrie_ilRL::level_pos = new uint64_t[height];

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
                    binTrie_ilRL::level_pos[level] = aux_level_pos[level];
                }
                else {
                    global_level_pos += aux_level_pos[level];
                    binTrie_ilRL::level_pos[level] = global_level_pos;
                }

            }

            binTrie_ilRL::bTrie = new bit_vector_il<block_size>(*_bTrie);
            delete _bTrie; 
            binTrie_ilRL::lastLevel = new bit_vector_il<block_size>(*_lastLevel);
            delete _lastLevel;
            binTrie_ilRL::b_rank = rank_support_il<1, block_size>(binTrie_ilRL::bTrie);

            binTrie_ilRL::lL_rank = rank_support_il<1, block_size>(binTrie_ilRL::lastLevel);
            
        };


        inline vector<uint64_t> * encodeRuns(std::vector<uint64_t> * posOnes) {
            vector<uint64_t> ones_to_write[binTrie_ilRL::height];
            vector<uint64_t> * ones_del_l_level = new vector<uint64_t>();
            uint64_t* aux_level_pos = new uint64_t[binTrie_ilRL::height];
            for (uint64_t i = 0; i < binTrie_ilRL::height; ++i) aux_level_pos[i] = 0;

            bool itsOneOne = false;
            binTrie_ilRL::writeCompressTrie(ones_to_write, ones_del_l_level, aux_level_pos, posOnes, 0, 0, itsOneOne);
            binTrie_ilRL::writeOnes(ones_to_write, aux_level_pos);

            delete[] aux_level_pos;   
            return ones_del_l_level;
        };

        
        inline void recursiveDecode(vector<uint64_t> &decoded, uint64_t partial_int, uint64_t node_id, uint16_t curr_level) {
            
            if (curr_level == binTrie_ilRL::height) {
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
            if (curr_level == binTrie_ilRL::height) {
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
            if (binTrie_ilRL::runs_encoded) {
                if (binTrie_ilRL::empty_trie) {
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
    
};

#endif