#include <iostream>
#include <vector>
#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/rank_support_v5.hpp>
#include "../src/intersection.hpp"
#include <sdsl/vectors.hpp>
#include <sdsl/int_vector.hpp>

#include "../src/flatBinTrieRL.hpp"
#include "../src/hybrid_structure_flat.hpp"
#include "../src/binTrie_ilRL.hpp"
#include "../src/hybrid_structure_il.hpp"

using namespace std;
using namespace sdsl;



bool runs = false;
bool verbose = false;
bool parallel = true;
uint32_t block_size = 512; //Only needed on binTrie_il
uint32_t low_p = 2;

template <class HybridStructure>
map<uint64_t, HybridStructure*> loadSequencesHS(std::ifstream &in, vector<vector<uint64_t>> &queries, u_int64_t n){
    vector<uint64_t> setIndexes;
    for(auto q: queries)
        setIndexes.insert(setIndexes.end(), q.begin(), q.end());
    sort(setIndexes.begin(), setIndexes.end());
    setIndexes.erase(unique(setIndexes.begin(), setIndexes.end()), setIndexes.end());
    map<uint64_t, HybridStructure*> sequences;

    uint64_t nIl = 0;
    uint32_t np = 0;
    for(uint32_t i = 0; i < n; ++i) {
        HybridStructure* trie = new HybridStructure();
        trie -> load(in);
        if (i == setIndexes[np]){
            sequences.insert({i, trie});
            np++;
            if (np == setIndexes.size()) break;
        }
        else { 
            trie -> free();
            delete trie;
        }
    }
    return sequences;
}

vector<vector<uint64_t>> loadQueryLog(std::string path, uint64_t numSets) {
    vector<vector<uint64_t>> queries;
    std::ifstream in(path);
    if (!in.is_open()) {
        cout << "Cant't open file: " << path << endl;
        return queries;
    }

    std::string line;
    while(std::getline(in, line)) {
        std::vector<uint64_t> query;
        std::istringstream iss(line);
        for (std::string s; iss >> s;) {
            uint64_t id = (uint64_t)stoull(s);
            if(id < numSets)
                query.push_back(id);
        }
        if(query.size() > 1)
            queries.push_back(query);
    }
    in.close();
    return queries;
}

template <class HybridStructure, uint8_t low_part_size>
void performIntersectionsHS( std::ifstream &in_sequences, std::string query_path,
                           std::string out_path, bool runs_encoded, bool intersectBF, bool parallel) {

    uint64_t trep = 10;
    vector<vector<uint64_t>> queries;
    map<uint64_t, HybridStructure*> sequences;

    uint32_t _1, u, n;
    in_sequences.read(reinterpret_cast<char*>(&n), sizeof(n));
    in_sequences.read(reinterpret_cast<char*>(&_1), sizeof(_1));
    in_sequences.read(reinterpret_cast<char*>(&u), sizeof(u));
    std::cout << "Num. of sets: " << n << std::endl;
    std::cout << "Universe: "<< u << std::endl;


    queries = loadQueryLog(query_path, n);
    cout << "Queries loaded succefully (" << queries.size()<<")" << endl;

    sequences = loadSequencesHS<HybridStructure>(in_sequences, queries, n);
    cout << "Sequences loaded succefully (" << sequences.size()<<")"<< endl;
    


    std::ofstream out;
    if (out_path != "") {
        out.open(out_path, std::ios::out);
        out << "elements_per_query,time execution,size_intersection" << std::endl; 
    }

    uint64_t nq = 0;
    long unsigned int total_time = 0;

    for (auto q: queries) {
        vector<HybridStructure> Hs;
        for (auto i: q){
            Hs.push_back(*sequences[i]);            
        }

        vector<uint64_t> intersection;
        if (Hs.size() <= 16){
            uint64_t time_10 = 0;
            // if(nq <= 11) {
                // trep = 1;
                for(int rep = 0; rep < trep; ++rep) {
                    auto start = std::chrono::high_resolution_clock::now();
                    Intersect_HS<HybridStructure, low_part_size>(Hs, intersection, intersectBF, runs_encoded, parallel);
                    auto end = std::chrono::high_resolution_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    auto time = elapsed.count();
                    total_time += time;
                    time_10 += time;
                    if (rep != trep-1)
                        intersection.clear();    
                }
                if (out.is_open()) {
                    out << Hs.size() << "," << (double)(time_10*1e-3)/trep<< "," << intersection.size()<< std::endl;
                }
            // }
            ++nq;
            if (nq % 1000 == 0 && verbose) {
                std::cout << nq << " queries processed" << std::endl;
            }
        }
    }

    std::cout << "Number of queries: " << nq << std::endl;
    std::cout <<"Avg time execution: " << (double)(total_time*1e-3)/(nq*trep) << "[ms]" << std::endl;

    out.close();
}
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il< select_support_scan<1>, 512, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il< select_support_mcl<1>, 512, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_scan<1>, 512, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);


template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il< select_support_mcl<1>, 512, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_scan<1>, 512, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_scan<1>, 512, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_scan<1>, 512, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS<HybridStructure_il<select_support_scan<1>, 512, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);


int main(int argc, char const *argv[]) {
    int mandatory = 5;
    // (*) mandatory parameters
    if (argc < mandatory){
        std::cout   << "collection filename "   // (*)
                        "query log "            // (*)
                        "[--parallel p] "       
                        "[--type_intersection type] " // (*)
                        "[--out r]"
                    <<
        std::endl;
        return 1;
    }


    int rank = 0;
    int select = 0;

    bool intersectBF = false;
    int type_intersection = 2;


    std::string sequences_filename = std::string(argv[1]);
    std::string querylog_filename   = std::string(argv[2]);
    std::string output_filename = "";
    
    ifstream in_sequences;
    in_sequences.open(sequences_filename, std::ios::binary | std::ios::in);

    for (int i = 1; i < argc; ++i){
        
        if (std::string(argv[i]) == "--parallel") {
            ++i;
            if (std::string(argv[i]) == "t")
                parallel = true;
            else
                parallel = false;
        }
        
        if (std::string(argv[i]) == "--type_intersection") {
            ++i;
            if (std::string(argv[i]) == "bk"){
                type_intersection = 1;
            }
            else if (std::string(argv[i]) == "bf"){
                type_intersection = 2;
                intersectBF = true;
            }
            else {
                type_intersection = 0;
            }
        }
        if (std::string(argv[i]) == "--out") {
            ++i;
            output_filename = std::string(argv[i]);
        }
        if (std::string(argv[i]) == "--verbose") {
            verbose = true;
        }
    }
    
    
    in_sequences.read(reinterpret_cast<char *>(&rank) ,sizeof(rank));
    if (rank == 1)
        in_sequences.read(reinterpret_cast<char *>(&block_size), sizeof(block_size));
    in_sequences.read(reinterpret_cast<char *>(&select), sizeof(select));
    in_sequences.read(reinterpret_cast<char *>(&runs), sizeof(runs));
    in_sequences.read(reinterpret_cast<char *>(&low_p), sizeof(low_p));


    std::cout << "Rank: ";
    if (rank == 0) std::cout << "rank v" << std::endl;
    else if (rank == 1) {
        std::cout << "rank il" << std::endl;
        std::cout << "Block size: " << block_size << std::endl;
    }
    else std::cout << "rank v5" << std::endl;
    std::cout << "Parallel:" << (parallel == 1 ? "true" : "false") << std::endl;
    std::cout << "Runs:" << (runs == 1 ? "true" : "false") << std::endl;
    std::cout << "Select: ";
    if (select  == 0){
        std::cout << "select mcl" << std::endl;
    }
    else if (select == 1) {
        std::cout << "select scan" << std::endl;
    } else std::cout << "-" << std::endl;
    std::cout << "Low part size: " << low_p << std::endl;

    std::cout << "Intersection algorithm: ";
    if (type_intersection == 2){
        std::cout << " Brute Force" << std::endl;
    }
    else if (type_intersection == 1) {
        std::cout << "Barbay and Kenyon" << std::endl;
    } 
    else
        std::cout << "-" << std::endl;

    std::cout << "outfilename: " << output_filename << endl;



    if(low_p == 1) {

            if (rank == 0 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } else if (rank == 1 && select == 0 ) {
                performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS<HybridStructure_il< select_support_scan<1>, 512, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
        
        
    } 
    else if(low_p == 2) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 0 ) {
                performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS<HybridStructure_il< select_support_scan<1>, 512, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
    }
    else if(low_p == 4) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } else if (rank == 1 && select == 0 ) {
                performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS<HybridStructure_il< select_support_scan<1>, 512, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
    } 
    else if(low_p == 8) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } else if (rank == 1 && select == 0 ) {
                performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS<HybridStructure_il< select_support_scan<1>, 512, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
        
    } 
    else if(low_p == 12) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 0 ) {
                performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS<HybridStructure_il< select_support_scan<1>, 512, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
        
    } 
    else if(low_p == 16) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 0 ) {
                performIntersectionsHS<HybridStructure_il<select_support_mcl<1>, 512, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS<HybridStructure_il< select_support_scan<1>, 512, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
        }
        
    
    
    std::cout << "END" << endl;

    return 0;
}