#include <iostream>
#include <vector>
#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/rank_support_v5.hpp>
// #include "../src/intersection.hpp"
#include "../src/flatBinTrie.hpp"
#include "../src/delta_trie_certificate.hpp"
#include <sdsl/vectors.hpp>
#include <sdsl/int_vector.hpp>



using namespace std;
using namespace sdsl;



template <class trieType>
vector<trieType>* loadSequences(std::string input_path){

    vector<trieType>* sequences = new vector<trieType>();
    std::ifstream in;
    in.open(input_path, std::ios::binary | std::ios::in);
    if (!in.is_open()) {
        cout << "Cant't open file: " << input_path << endl;
        return sequences;
    }

    uint32_t _1, u, n;
    // Read universe
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    in.read(reinterpret_cast<char*>(&_1), sizeof(_1));
    in.read(reinterpret_cast<char*>(&u), sizeof(u));
    std::cout << "Num. of sets: " << n << std::endl;
    std::cout << "Universe: "<< u << std::endl;
    
    uint64_t count = 0;
    while ( count < n) { 
        if (in.eof()){
            break;
        }
        trieType b = trieType();
        b.load(in);
        sequences -> push_back(b);
        count++;
    }
    return sequences;
}
template 
vector<flatBinTrie<sdsl::rank_support_v<1>>>* loadSequences<flatBinTrie<sdsl::rank_support_v<1>>>(std::string input_path);
template 
vector<flatBinTrie<sdsl::rank_support_v5<1>>>* loadSequences<flatBinTrie<sdsl::rank_support_v5<1>>>(std::string input_path);



vector<vector<uint32_t>>* loadQueryLog(std::string path, uint32_t numSets) {
    vector<vector<uint32_t>>* queries = new vector<vector<uint32_t>>();
    std::ifstream in(path);
    if (!in.is_open()) {
        cout << "Cant't open file: " << path << endl;
        return queries;
    }

    std::string line;
    while(std::getline(in, line)) {
        std::vector<uint32_t> query;

        std::istringstream iss(line);
        for (std::string s; iss >> s;) {
            uint32_t id = (uint32_t)stoi(s);
            if(id < numSets)
                query.push_back(id);
        }
        if(query.size() > 1)
            queries->push_back(query);
    }
    
    return queries;
}

template <class trieType>
void performDeltaTrie( std::string sequences_path, std::string query_path,
                           std::string out_path, bool runs_encoded) {

    uint64_t trep = 1;
    std::ofstream out;
    if (out_path != "") {
        out.open(out_path, std::ios::out);
        out << "elements_per_query,time execution,size_intersection,query number, delta" << std::endl; 
    }

    vector<vector<uint32_t>>* queries;
    vector<trieType>* sequences;

    sequences = loadSequences<trieType>(sequences_path);
    cout << "Sequences loaded succefully (" << sequences->size()<<")"<< endl;
    queries   = loadQueryLog(query_path, sequences->size());
    cout << "Queries loaded succefully (" << queries->size()<<")" << endl;

    uint64_t n_ranks = 0;

    uint64_t nq = 0;
    long unsigned int total_time = 0;
    long unsigned int total_delta = 0;
    uint64_t total_height = 0;
    uint64_t time_of_ranks = 0;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < queries -> size(); ++i) {
    // for (int i = 2; i < 3; ++i) {
    //     std::cout << endl << "queries: " << i << std::endl;
        vector<trieType> Bs;
        Bs.reserve(16);
        for(int j = 0; j < (*queries)[i].size(); ++j) {
            Bs.push_back((*sequences)[(*queries)[i][j]]);
        }
        vector<uint64_t> intersection;
        uint64_t delta = 0;
        if (Bs.size() <= 16){
            uint64_t time_10 = 0;
            for(int rep = 0; rep < trep; ++rep) {
                delta = 0;
                auto start = std::chrono::high_resolution_clock::now();
                DeltaTrieCertificate(Bs, intersection,delta, runs_encoded);
                auto end = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                auto time = elapsed.count();
                total_time += time;
                total_delta += delta;
                time_10 += time;
                if (rep != trep-1)
                    intersection.clear();        
            }
            if (out.is_open()) {
                out << Bs.size() << "," << (float)time_10/trep<< "," << intersection.size()<< "," <<i<< "," <<delta<< std::endl;
            }
            // std::cout << nq <<"|Time execution: " << (double)(time_10*1e-3)/(trep) << "[ms] " << intersection.size() << std::endl;
            ++nq;
            if (nq % 1000 == 0) {
                std::cout << nq << " queries processed" << std::endl;
            }
        }
    }

    std::cout << "Number of queries: " << nq << std::endl;
    std::cout <<"Avg time execution: " << (double)(total_time*1e-3)/(nq*trep) << "[ms]" << std::endl;
    std::cout <<"Avg delta trie certificate: " << (double)(total_delta)/(nq*trep) << std::endl;

    out.close();

    delete sequences;
    delete queries;
}

template
void performDeltaTrie<flatBinTrie<sdsl::rank_support_v<1>>>(std::string sequences_path, std::string query_path, std::string out_path,bool runs_encoded);
template
void performDeltaTrie<flatBinTrie<sdsl::rank_support_v5<1>>>(std::string sequences_path, std::string query_path, std::string out_path,bool runs_encoded);




int main(int argc, char const *argv[]) {
    int mandatory = 5;
    // (*) mandatory parameters
    if (argc < mandatory){
        std::cout   << "collection filename "   // (*)
                        "query log "            // (*)
                        "[--runs r] "           // (*)
                        "[--rank v] "
                        "[--out results_filename]"
                    <<
        std::endl;
        return 1;
    }

    int rank = 0;
    int select = 0;
    uint32_t block_size = 512;
    bool runs = true;
    std::string sequences_filename = std::string(argv[1]);
    std::string querylog_filename   = std::string(argv[2]);
    std::string output_filename = "";
    for (int i = 1; i < argc; ++i){
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
    }
    
    std::cout << "Rank: ";
    if (rank == 0) std::cout << "rank v" << std::endl;
    else std::cout << "rank v5" << std::endl;
    std::cout << "Runs:" << (runs == 1 ? "true" : "false") << std::endl;
    

    std::cout << "outfilename: " << output_filename << endl;

        if (rank == 0) {
            performDeltaTrie<flatBinTrie<sdsl::rank_support_v<1>>>(sequences_filename, querylog_filename, output_filename, runs);
        } 
        else {
            performDeltaTrie<flatBinTrie<sdsl::rank_support_v5<1>>>(sequences_filename, querylog_filename, output_filename, runs);
        }

    
    std::cout << "END" << endl;

    return 0;
}
