#include <iostream>
#include <fstream>
#include "../src/flatBinTrieRL.hpp"
#include "../src/hybrid_structure_flat.hpp"
#include <sdsl/bit_vectors.hpp>
#include "../src/binTrie_ilRL.hpp"
#include "../src/hybrid_structure_il.hpp"



using namespace std;
using namespace sdsl;

bool verbose = false;
uint32_t block = 512;
uint8_t low_p = 2;

vector<uint64_t>* read_inv_list(std::ifstream &input_stream, uint32_t n) {

    vector <uint64_t>* il = new vector<uint64_t>();
    il -> reserve(n);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t x;
        input_stream.read(reinterpret_cast<char *>(&x), 4);
        il -> push_back((uint64_t)x);
    }
    return il;
}

template <uint32_t block_size, uint8_t low_part_size>
void buildCollection(std::string input_path, std::string out_path, uint64_t min_size,
                    uint64_t max_size, int rank_type, bool runs, int select_type) {
    std::ifstream input_stream;
    input_stream.open(input_path, std::ios::binary | std::ios::in);
    if (!input_stream.is_open()) {
        cout << "Can't open input file:  " << input_path << endl;
        return;
    }
    std::ofstream out;
    if (out_path != ""){
        out.open(out_path, ios::out | ios::binary);
        if(!out.is_open()) {
            cout << "Can't open  output file:  " << out_path << endl;
            return;
        }
    }
    

    u_int64_t _1, u;
    // Read universe of collection
    input_stream.read(reinterpret_cast<char *>(&_1), sizeof(_1));
    input_stream.read(reinterpret_cast<char *>(&u), sizeof(u));

    u_int64_t nSets = 0;
    while (true) {
        uint32_t n;
        input_stream.read(reinterpret_cast<char *>(&n), 4);
        if (input_stream.eof()){
            break;
        }
        if (n > min_size && (n <= max_size || max_size == 0)) {
            nSets++;
        }
        input_stream.seekg(4*n, ios::cur);
    }
    // Set file pointer in first set
    input_stream.clear();
    input_stream.seekg(2*sizeof(u), ios::beg);
    // uint32_t low_p = low_part_size;
    // Write universe in out file
    if (out_path != "") {
        out.write(reinterpret_cast<char *> (&rank_type), sizeof(rank_type));
        if (rank_type == 1)
            out.write(reinterpret_cast<char *> (&block), sizeof(block));
        out.write(reinterpret_cast<char *> (&select_type), sizeof(select_type));
        out.write(reinterpret_cast<char *> (&runs), sizeof(runs));
        out.write(reinterpret_cast<char *> (&low_p), sizeof(low_p));
        out.write(reinterpret_cast<char *> (&nSets), sizeof(nSets));
        out.write(reinterpret_cast<char *> (&_1), sizeof(_1));
        out.write(reinterpret_cast<char *> (&u), sizeof(u));
    }
    cout << "Universe: " << u << endl;

    uint64_t total_size = 0;
    uint64_t total_hs_size = 0;
    uint64_t total_elements = 0;
    uint64_t n_il = 0;

    uint32_t max_n = 0;
    while (true) {        
        uint32_t n;
        input_stream.read(reinterpret_cast<char *>(&n), 4);
        if (input_stream.eof()){
            break;
        }
        
        uint64_t trie_bytes_size = 0;
        uint64_t hs_bytes_size = 0;
        if (n > min_size && (n <= max_size || max_size == 0)){
            vector <uint64_t> *il = read_inv_list(input_stream, n);
            uint64_t max_value = (*il)[n - 1];
            max_n = max_n < n ? n : max_n;
            
            if (rank_type == 0 && select_type == 0 ) {
                HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, low_part_size> hs;
                hs = HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, low_part_size>(*il, u);
                if (runs)
                    hs.encodeRuns();
                if (out_path != "") {
                    hs.serialize(out);
                }
                hs_bytes_size = hs.size_in_bytes();
                hs.free();
            }
            else if (rank_type == 0  && select_type == 1 ) {
                HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, low_part_size> hs;
                hs = HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, low_part_size>(*il, u);
                if (runs)
                    hs.encodeRuns();
                if (out_path != "") {
                    hs.serialize(out);
                }
                hs_bytes_size = hs.size_in_bytes();
                hs.free();                
            }
            else if (rank_type == 2 && select_type == 0 ) {
                HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, low_part_size> hs;
                hs = HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, low_part_size>(*il, u);
                if (runs)
                    hs.encodeRuns();
                if (out_path != "") {
                    hs.serialize(out);
                }
                hs_bytes_size = hs.size_in_bytes();
                hs.free();
            }
            else if (rank_type == 2  && select_type == 1 ) {
                HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, low_part_size> hs;
                hs = HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, low_part_size>(*il, u);
                if (runs)
                    hs.encodeRuns();
                if (out_path != "") {
                    hs.serialize(out);
                }
                hs_bytes_size = hs.size_in_bytes();
                hs.free();                
            }
            else if (rank_type == 1 && select_type == 0 ) {
                HybridStructure_il< select_support_mcl<1>, block_size, low_part_size> hs;
                hs = HybridStructure_il< select_support_mcl<1>, block_size, low_part_size>(*il, u);
                if (runs)
                    hs.encodeRuns();
                if (out_path != "") {
                    hs.serialize(out);
                }
                hs_bytes_size = hs.size_in_bytes();
                hs.free();
            }
            else if (rank_type == 1 && select_type == 1 ) {
                HybridStructure_il< select_support_scan<1>, block_size, low_part_size> hs;
                hs = HybridStructure_il<select_support_scan<1>, block_size, low_part_size>(*il, u);
                if (runs)
                    hs.encodeRuns();
                if (out_path != "") {
                    hs.serialize(out);
                }
                hs_bytes_size = hs.size_in_bytes();
                hs.free();
            }            

            total_size += trie_bytes_size;
            total_elements += n;
            n_il++;
            total_hs_size += hs_bytes_size;

            if ((n_il % 1000) == 0 && verbose) {
                cout << n_il  <<" Sets processed " << endl;
            }

            delete il;
        }

        else {
            input_stream.seekg(4*n, ios::cur);
        }
        
    }
    input_stream.close();
    out.close();

    cout << "Maximum number of elements in sets: " << max_n << endl;
    cout << "Total inverted lists: " << n_il << "| Bpi: " << (float)(total_size*8)/total_elements << endl;
    cout << "Total ints: " << total_elements << endl;
    cout << "Total HS: " << n_il << "| Bpi: " << (float)(total_hs_size*8)/total_elements << endl;
    
    return;
}


int main(int argc, char** argv) {
    int mandatory = 3;

    if (argc < mandatory){
        std::cout   << "collection filename (*)"       // (*)
                        "[--min_size m] "           // (*)
                        "[--max_size m] "           // (*)
                        "[--rank v] (*)"               // (*)
                        "[--runs r] (*)"               // (*)
                        "[--low_part_size l] (*)"      // (*)
                        "[--out output_filename]"
                    <<
        std::endl;
        return 1;
    }

    int rank = 0;
    int select = 0;
    uint64_t min_size = 0;
    uint64_t max_size = 0;
    // uint32_t block_size = 512;
    // uint32_t low_part_size = 2;
    
    bool runs = false;
    std::string output_filename = "";
    std::string input_filename = std::string(argv[1]);
    std::cout << input_filename << endl;
    for (int i = 2; i < argc; ++i){
        if (std::string(argv[i]) == "--min_size") {
            ++i;
            min_size = std::stoull(argv[i]);
        }
        if (std::string(argv[i]) == "--max_size") {
            ++i;
            max_size = std::stoull(argv[i]);
        }
        if (std::string(argv[i]) == "--rank") {
            ++i;
            if (std::string(argv[i]) == "v") {
                rank = 0;
            }
            else if (std::string(argv[i]) == "il") {
                rank = 1;
                i++;
                block = std::atoi(argv[i]);
                // block = block != 128 && block != 256 ? 512 : block;
            }
            else {
                rank = 2;
            }
        }
        if (std::string(argv[i]) == "--runs") {
            ++i;
            if (std::string(argv[i]) == "t")
                runs = true;
            else
                runs = false;
        }
        // if (std::string(argv[i]) == "--select") {
        //     ++i;
        //     if (std::string(argv[i]) == "mcl")
        //         select = 0;
        //     else if (std::string(argv[i]) == "scan")
        //         select = 1;
        //     else 
        //         select = 2;    
        // }
        if (std::string(argv[i]) == "--low_part_size") {
            ++i;
            low_p = std::atoi(argv[i]);
        }
        if (std::string(argv[i]) == "--verbose") {
            ++i;
            if (std::string(argv[i]) == "t")
                verbose = true;
            else
                verbose = false;
        }
        if (std::string(argv[i]) == "--out") {
            ++i;
            output_filename = std::string(argv[i]);
        }
    }
    
    std::cout << "Min size: " << min_size << std::endl;
    std::cout << "Max size: " << max_size << std::endl;
    std::cout << "Rank: ";
    if (rank  == 0){
        std::cout << "rank v" << std::endl;
    }
    else if (rank == 1) {
        std::cout << "rank il" << std::endl;
        std::cout << "Block size:" << block << std::endl;
    }
    else {
        std::cout << "rank v5" << std::endl;
    }
    std::cout << "Runs: " << (runs == 1 ? "true" : "false") << std::endl;
   
    std::cout << "Low part size: " << low_p << std::endl;
    
    std::cout << "Output file name: "<< output_filename << endl;
    
    switch (low_p)  {
        case 1:
            if (block == 128)
                buildCollection<128, 1>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else if(block == 256)
                buildCollection<256, 1>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else 
                buildCollection<512, 1>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            break;
        case 2:
            if (block == 128)
                buildCollection<128, 2>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else if(block == 256)
                buildCollection<256, 2>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else   
                buildCollection<512, 2>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            break;
        case 4:
            if (block == 128)
                buildCollection<128, 4>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else if(block == 256)
                buildCollection<256, 4>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else   
                buildCollection<512, 4>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            break;
        case 8:
            if (block == 128)
                buildCollection<128, 8>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else if(block == 256)
                buildCollection<256, 8>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else    
                buildCollection<512, 8>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            break;
        case 12:
            if (block == 128)
                buildCollection<128, 12>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else if(block == 256)
                buildCollection<256, 12>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else
                buildCollection<512, 12>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            break;
        case 16:
        default:
            if (block == 128)
                buildCollection<128, 16>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else if(block == 256)
                buildCollection<256, 16>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            else
                buildCollection<512, 16>(input_filename, output_filename, min_size, max_size, rank, runs, select);
            break;
    }
}