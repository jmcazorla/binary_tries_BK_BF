#include <iostream>
#include <fstream>
#include "../src/flatBinTrie.hpp"
#include <sdsl/bit_vectors.hpp>


using namespace std;
using namespace sdsl;


bool set_num_64 = false;

vector<uint64_t>* read_inv_list_64(std::ifstream &input_stream, uint32_t n) {
    vector <uint64_t>* il = new vector<uint64_t>();
    il -> reserve(n);
    for (uint32_t i = 0; i < n; ++i) {
        uint64_t x;
        input_stream.read(reinterpret_cast<char *>(&x), sizeof(uint64_t));
        il -> push_back(x);
    }
    return il;
}

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

template <uint32_t block_size>
void buildCollection(std::string input_path, std::string out_path, uint64_t min_size, uint64_t max_size, int rank_type, bool runs) {
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
    

    uint32_t _1, u;
    // Read universe of collection
    input_stream.read(reinterpret_cast<char *>(&_1), sizeof(_1));
    // input_stream.read(reinterpret_cast<char *>(&u), sizeof(u));
    set_num_64 ? input_stream.read(reinterpret_cast<char *>(&u), sizeof(uint64_t)) : input_stream.read(reinterpret_cast<char *>(&u), sizeof(uint32_t));

    uint32_t nSets = 0;
    while (true) {
        uint32_t n;
        input_stream.read(reinterpret_cast<char *>(&n), 4);
        if (input_stream.eof()){
            break;
        }
        if (n > min_size && (n <= max_size || max_size == 0)) {
            nSets++;
        }

        input_stream.seekg(n*(set_num_64 ? sizeof(uint64_t) : sizeof(uint32_t)), ios::cur);
    }
    // Set file pointer in first set
    input_stream.clear();

    set_num_64 ? input_stream.seekg(sizeof(uint32_t) + sizeof(uint64_t), ios::beg) : input_stream.seekg(2*sizeof(uint32_t), ios::beg);
    // Write universe in out file
    if (out_path != "") {
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
    uint32_t print = 1;
    
    while (true) {
        
        uint32_t n;
        input_stream.read(reinterpret_cast<char *>(&n), 4);
        if (input_stream.eof()){
            break;
        }
        
        uint64_t trie_bytes_size = 0;
        if (n > min_size && (n <= max_size || max_size == 0)){
            vector <uint64_t> * il = set_num_64 ?  read_inv_list_64(input_stream, n) : read_inv_list(input_stream, n);
            uint64_t max_value = (*il)[n - 1];
            max_n = max_n < n ? n : max_n;

            if(print == 0) cout  <<rank_type<< " - "<< " - " << endl;


            if (rank_type == 0 ) {

                    flatBinTrie<rank_support_v<1>> trie_v;
                    trie_v = flatBinTrie<rank_support_v<1>>(*il, u);
                    if (runs)
                        trie_v.encodeRuns();
                    if (out_path != "") {
                        trie_v.serialize(out);
                    }
                    trie_bytes_size = trie_v.size_in_bytes();
                    trie_v.free(); 
               
            }
            if (rank_type == 1  ) {

                    flatBinTrie<rank_support_v5<1>> trie_v5;
                    trie_v5 = flatBinTrie<rank_support_v5<1>>(*il, u);
                    if (runs)
                        trie_v5.encodeRuns();
                    if (out_path != "") {
                        trie_v5.serialize(out);
                    }
                    trie_bytes_size = trie_v5.size_in_bytes();
                    trie_v5.free(); 
                
            }


            total_size += trie_bytes_size;
            total_elements += n;
            n_il++;

            // cout << "#Elements: " << n << " | Bpi: " << (float)(trie_bytes_size*8)/n << endl;
            if ((n_il % 1000) == 0) {
                cout << n_il  <<" Sets processed " << endl;
            }
            delete il;

        }
        else {
            input_stream.seekg(n*(set_num_64 ? sizeof(uint64_t) : sizeof(uint32_t)), ios::cur);
        }
       
    }
    input_stream.close();
    out.close();


    cout << "Maximum number of elements in sets: " << max_n << endl;
    cout << "Total inverted lists: " << n_il << "| Bpi: " << (float)(total_size*8)/total_elements << endl;
    cout << "Total ints: " << total_elements << endl;
    
    return;
}


int main(int argc, char** argv) {
    int mandatory = 3;

    if (argc < mandatory){
        std::cout   << "collection filename "       // (*)
                        "[--min_size m] "           // (*)
                        "[--max_size m] "           // (*)
                        "[--rank v] "               // (*)
                        "[--runs r] "               // (*)
                        "[--out output_filename]"
                    <<
        std::endl;
        return 1;
    }

    int rank = 0;
    uint64_t min_size = 0;
    uint64_t max_size = 0;
    uint32_t block_size = 512;

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
            else {
                rank = 1;
            }
        }
        if (std::string(argv[i]) == "--runs") {
            ++i;
            if (std::string(argv[i]) == "t")
                runs = true;
            else
                runs = false;
        }
        if (std::string(argv[i]) == "--out") {
            ++i;
            output_filename = std::string(argv[i]);
        }
        if (std::string(argv[i]) == "--set_num_64")
            set_num_64 = true;
    }
    
    std::cout << "Min size: " << min_size << std::endl;
    std::cout << "Max size: " << max_size << std::endl;
    std::cout << "Rank: ";
    if (rank  == 0){
        std::cout << "rank v" << std::endl;
    }
    else {
        std::cout << "rank v5" << std::endl;
    }
    std::cout << "Runs: " << (runs == 1 ? "true" : "false") << std::endl;
    
    std::cout << " set_num_64: " << (set_num_64 == 1 ? "true" : "false") << std::endl;
    std::cout << "Output file name: "<< output_filename << endl;
    
    // Call function here
    
    if (block_size == 512) {
            buildCollection<512>(input_filename, output_filename, min_size, max_size, rank, runs);
    } else if (block_size == 256) {
            buildCollection<256>(input_filename, output_filename, min_size, max_size, rank, runs);
    } else {
            buildCollection<512>(input_filename, output_filename, min_size, max_size, rank, runs);
    }

        
    // std::cout << "Presiona Enter para sainput_streamlir...";
    // std::cin.get(); // Espera a que el usuario presione Enter
    // return 0;
}