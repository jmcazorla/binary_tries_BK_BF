#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/rank_support_v5.hpp>
#include "../src/fastBinaryTrie.hpp"
#include "../src/new_baseline.hpp"
#include "../src/barbay_and_kenyon.hpp"

using namespace std;
using namespace sdsl;

bool runs = false;
bool levelwise = false;
uint32_t block_size = 512; //Only needed on binTrie_il
uint16_t wsize = 64;
bool verbose = false;

enum class AlgoType {
    SWAPPING_SVS,
    BAEZA_YATES,
    SMALL_ADAPTIVE,
    BARBAY_KENYON
};

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

template<class rankType, class wordType>
map<uint64_t, fastBinaryTrie<rankType, wordType>*> loadTries(ifstream &in, vector<vector<uint64_t>> &queries, uint64_t numSets){
    vector<uint64_t> setIndexes;
    for(auto q: queries)
        setIndexes.insert(setIndexes.end(), q.begin(), q.end());
    sort(setIndexes.begin(), setIndexes.end());
    setIndexes.erase(unique(setIndexes.begin(), setIndexes.end()), setIndexes.end());

    map<uint64_t, fastBinaryTrie<rankType, wordType>*> tries; 
    uint64_t nIl = 0;
    uint32_t np = 0;
    for(uint32_t i = 0; i < numSets; ++i) {
        fastBinaryTrie<rankType, wordType>* trie = new fastBinaryTrie<rankType, wordType>();
        trie -> load(in);
        if (i == setIndexes[np]){
            tries.insert({i, trie});
            np++;
            if (np == setIndexes.size()) break;
        }
        else { // I dont now how space will use the trie
            delete trie;
        }
    }
    return tries;
}

template<class rankType, class wordType>
void performIntersections( std::ifstream &in_sequences, std::string query_path,
                         std::string out_path, bool runs_encoded, AlgoType algo) {
    
    uint64_t trep = 10;
    vector<vector<uint64_t>> queries;
    map<uint64_t, fastBinaryTrie<rankType, wordType>*> tries;

    uint32_t _1, u, n;
    in_sequences.read(reinterpret_cast<char*>(&n), sizeof(n));
    in_sequences.read(reinterpret_cast<char*>(&_1), sizeof(_1));
    in_sequences.read(reinterpret_cast<char*>(&u), sizeof(u));
    std::cout << "Num. of sets: " << n << std::endl;
    std::cout << "Universe: "<< u << std::endl;

    queries = loadQueryLog(query_path, n);
    cout << "Queries loaded succefully, Total: " << queries.size() << "" << endl;
    tries = loadTries<rankType, wordType>(in_sequences, queries, n);
    cout << "Sequences loaded succefully, Total: " << tries.size() << endl;

    std::ofstream out;
    if (out_path != "") {
        out.open(out_path, std::ios::out);
        out << "elements_per_query,time execution,size_intersection" << std::endl; 
    }

    cout << "Computing queries...\n";
    uint64_t nq = 0;
    long unsigned int total_time = 0;
    for (auto q: queries) {
        vector<fastBinaryTrie<rankType, wordType>*> QTries;
        vector<vector<uint64_t>> Qsets;
        for (auto i: q){
            vector<uint64_t> decodedTrie;
            tries[i]->decode(decodedTrie);
            Qsets.push_back(decodedTrie);
            QTries.push_back(tries[i]);    
        }
        vector<uint64_t> intersection;
        if (QTries.size() <= 16){
            uint64_t time_10 = 0;
            // trep = 1;
            // if(nq <= 11) {
                trep = 1;
                for(int rep = 0; rep < trep; ++rep) {
                    auto start = std::chrono::high_resolution_clock::now();
                    // intersection = swapping_svs(Qsets);
                    // intersection = so_baeza_yates(Qsets);
                    // intersection = small_adaptive(Qsets);                    
                    // barbayKenyon(Qsets, Qsets.size(), intersection);
                    switch(algo) {
                        case AlgoType::SWAPPING_SVS:
                            intersection = swapping_svs(Qsets);
                            break;
                        case AlgoType::BAEZA_YATES:
                            intersection = so_baeza_yates(Qsets);
                            break;
                        case AlgoType::SMALL_ADAPTIVE:
                            intersection = small_adaptive(Qsets);
                            break;
                        case AlgoType::BARBAY_KENYON:
                            barbayKenyon(Qsets, Qsets.size(), intersection);
                            break;
                    }

                    auto end = std::chrono::high_resolution_clock::now();
                    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                    auto time = elapsed.count();
                    total_time += time;  
                    time_10 += time;  
                }
                if (out.is_open()) {
                    out << Qsets.size() << "," << (double)(time_10*1e-3)/trep<< "," << intersection.size()<< std::endl;
                }
            // }

            ++nq;
            if (nq % 1000 == 0 && verbose) {
                std::cout << nq << " queries processed" << std::endl;
            }
        }
    }
    for (auto T: tries)
        delete T.second;

    std::cout << "Number of queries: " << nq << std::endl;
    std::cout <<"Avg time execution: " << (double)(total_time*1e-3)/(nq) << "[ms]" << std::endl;
}


int main(int argc, char const *argv[]) {
    int mandatory = 3;
    AlgoType selected_algo = AlgoType::BARBAY_KENYON; // Default
    // (*) mandatory parameters
    if (argc < mandatory){
        std::cout   << "collection filename " // (*)
                        "query log" // (*)
                        "[--algo a]" //(svs, by, adaptive, bk)
                        "[--out r]"
                    << "\n";
        return 1;
    }

    int rankT = 0;
    uint64_t min_size;
    std::string sequences_filename = std::string(argv[1]);
    std::string querylog_filename   = std::string(argv[2]);
    std::string output_filename = "";

    for (int i = 1; i < argc; ++i){

        if (std::string(argv[i]) == "--out") {
            ++i;
            output_filename = std::string(argv[i]);
        }
        if (std::string(argv[i]) == "--algo") {
            ++i;
            std::string algo = std::string(argv[i]);
            if (algo == "svs") selected_algo = AlgoType::SWAPPING_SVS;
            else if (algo == "by") selected_algo = AlgoType::BAEZA_YATES;
            else if (algo == "adaptive") selected_algo = AlgoType::SMALL_ADAPTIVE;
            else if (algo == "bk") selected_algo = AlgoType::BARBAY_KENYON;
        }
        if (std::string(argv[i]) == "--verbose") {
            verbose = true;
        }
    }

    ifstream in_sequences;
    in_sequences.open(sequences_filename, std::ios::binary | std::ios::in);

    in_sequences.read(reinterpret_cast<char *>(&rankT) ,sizeof(rankT));
    if (rankT == 1) 
        in_sequences.read(reinterpret_cast<char *>(&block_size), sizeof(block_size));
    in_sequences.read(reinterpret_cast<char *>(&runs), sizeof(runs));
    in_sequences.read(reinterpret_cast<char *>(&levelwise), sizeof(levelwise));
    in_sequences.read(reinterpret_cast<char *>(&wsize), sizeof(wsize));



    std::cout << "Type of trie\n";
    std::cout << "* Rank: ";
    if (rankT == 0) std::cout << "rank v" << std::endl;
    else if (rankT == 1) {
        std::cout << "rank il" << std::endl;
        std::cout << "** Block size: " << block_size << std::endl;
    }
    else std::cout << "rank v5" << std::endl;
    std::cout << "* Runs:" << (runs == 1 ? "true" : "false") << std::endl;
    std::cout << "* Level-wise: " << (levelwise == 1 ? "true" : "false") << std::endl;

    if (rankT == 0){
        if (wsize == 64)
            performIntersections<sdsl::rank_support_v<1>, uint64_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
        else if (wsize == 32)
            performIntersections<sdsl::rank_support_v<1>, uint32_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
        else if (wsize == 16)
            performIntersections<sdsl::rank_support_v<1>, uint16_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
        else
            performIntersections<sdsl::rank_support_v<1>, uint8_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
    }
    else{
        if (wsize == 64)
            performIntersections<sdsl::rank_support_v5<1>, uint64_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
        else if (wsize == 32)
            performIntersections<sdsl::rank_support_v5<1>, uint32_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
        else if (wsize == 16)
            performIntersections<sdsl::rank_support_v5<1>, uint16_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
        else
            performIntersections<sdsl::rank_support_v5<1>, uint8_t>(in_sequences, querylog_filename, output_filename, runs, selected_algo);
    }
        
    return 0;
}
