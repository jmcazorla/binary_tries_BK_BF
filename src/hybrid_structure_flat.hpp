#ifndef HYBRID_STRUCTURE_FLAT_HPP
#define HYBRID_STRUCTURE_FLAT_HPP

#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <sdsl/int_vector.hpp>
#include "flatBinTrieRL.hpp"
using namespace std;
using namespace sdsl;


template <class rankType = rank_support_v<1>, class selectType = select_support_mcl<1>, uint8_t low_part_size = 2>
class HybridStructure_flat {
private:
    flatBinTrieRL<rankType> *binTrie;
    sdsl::int_vector<low_part_size> *lowPart;
    sdsl::bit_vector *flags;
    selectType f_select;

public:
    HybridStructure_flat() = default;

    ~HybridStructure_flat() = default;

    HybridStructure_flat(std::vector<uint64_t>& set, uint64_t u) {
        uint64_t length_it = set.size();
        lowPart = new sdsl::int_vector<low_part_size>(length_it);
        flags = new sdsl::bit_vector(length_it);

        uint64_t num, highAux, lowAux;
        size_t j = 0;
        std::vector<uint64_t> highPartAux(length_it);
        for (size_t i = 0; i < length_it; i++) {
            (*flags)[i] = 0;
            num = set[i];
            highAux = num >> low_part_size;   // Top (shift right)
            lowAux = num & ((1 << low_part_size) - 1);   // Low part (bit mask)
            (*lowPart)[i] = lowAux;
            // Store upper part without repeating elements
            if (j == 0) {
                highPartAux[j++] = highAux;
            } else if (highPartAux[j - 1] != highAux) {
                highPartAux[j++] = highAux;
                (*flags)[i - 1] = 1;
            }
        }
        (*flags)[length_it - 1] = 1;
        std::vector<uint64_t> highPart(highPartAux.begin(), highPartAux.begin() + j);
        uint16_t height = floor(log2(u - 1)) +  1;
        uint64_t newU = (pow(2, (height - low_part_size))) - 1;
        HybridStructure_flat::binTrie = new flatBinTrieRL<rankType>(highPart, newU);
        HybridStructure_flat::f_select = selectType(flags);
    }
    
    enum { fixed_low_part_size = low_part_size }; 

    inline void free(){
        delete binTrie;
        delete lowPart;
        delete flags;

    }

    inline flatBinTrieRL<rankType> *getBinTrie() {
        return binTrie;
    }

    inline sdsl::int_vector<low_part_size> *getLowPart() {
        return lowPart;
    }

    inline sdsl::bit_vector *getFlags() {
        return flags;
    }

    inline uint64_t selectFlags(uint64_t value) {
        uint64_t select = HybridStructure_flat::f_select(value);
        return select;
    }

    inline uint64_t size_in_bytes() {
        uint64_t lowPart_size = sdsl::size_in_bytes(*lowPart);
        uint64_t flags_size = sdsl::size_in_bytes(*flags);
        uint64_t f_select_size = sdsl::size_in_bytes(f_select);
        return binTrie->size_in_bytes()
                        + lowPart_size
                        + flags_size
                        + f_select_size;
    }

    inline uint64_t serialize(std::ostream &out) {
        uint64_t bvs_size, select_size;
        bvs_size  = binTrie->serialize(out) + flags->serialize(out) + lowPart->serialize(out);
        select_size = f_select.serialize(out);

        return bvs_size + select_size;
    }
    
    inline void load(std::istream &in){
        binTrie = new flatBinTrieRL<rankType>();
        lowPart = new sdsl::int_vector<low_part_size>();
        flags = new sdsl::bit_vector();
        
        binTrie->load(in);
        flags->load(in);
        lowPart->load(in);
        f_select.load(in, flags);
    }

    inline void encodeRuns() {
        uint64_t chain_length = (1 << HybridStructure_flat::fixed_low_part_size);
        std::vector<uint64_t> * posOnes = new std::vector<uint64_t>(); // Vector to store the positions        
        uint64_t number_ones = HybridStructure_flat::binTrie->elements_coded();
        bool delete_prev = false;
        bool zero_prev = false;
        uint64_t pos_one = 1;

        uint64_t prev_one = HybridStructure_flat::f_select(1);
        if(prev_one == (chain_length-1)){
            posOnes->push_back(1);
            zero_prev = true;
            delete_prev = true;
        }
        for (uint64_t i = 2; i <=  number_ones; ++i){
            pos_one = HybridStructure_flat::f_select(i);
            if((prev_one+chain_length) == pos_one){
                posOnes->push_back(i);
                delete_prev = !zero_prev;
                zero_prev = true;
            } else {
                if(delete_prev) {
                    posOnes->pop_back();
                }
                zero_prev = false;
                delete_prev = false;
            }
            prev_one = pos_one;
        }
        if(delete_prev){
            posOnes->pop_back();
        }
        vector<uint64_t> * onesDelLLevel = HybridStructure_flat::binTrie->encodeRuns(posOnes);
        if( onesDelLLevel->size() > 0) {            
            uint64_t amount_delete = chain_length*onesDelLLevel->size();
            uint64_t length = lowPart->size() - amount_delete;
            uint64_t posDel = 0;
            long int posOk = 0;
            uint64_t ones = 0;

            posDel = 0;
            uint64_t iEnd = 0;
            posOk = -1;
            int zero = 0;
            for (iEnd = 0; iEnd < lowPart->size(); iEnd++)
            {
                posOk++;
                (*lowPart)[posOk] = (*lowPart)[iEnd];
                (*flags)[posOk] = (*flags)[iEnd];

                if((*flags)[iEnd] == 1) {
                    ones++;
                }

                if(ones == (*onesDelLLevel)[posDel] && posDel < onesDelLLevel->size()) {    
                    posDel++;                    
                    posOk = posOk - chain_length;
                }    

                if((*flags)[iEnd] == 1) {
                    zero = 0;
                } 
                else {
                    zero++;
                }                
            }
            posOk++;
            lowPart->resize(posOk);
            flags->resize(posOk);
            HybridStructure_flat::f_select = selectType(flags);            
        }
        delete onesDelLLevel;
        delete posOnes;
    }

};

#endif  // HYBRID_STRUCTURE_FLAT_HPP
