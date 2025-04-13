#include <iostream>
#include <vector>
#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/rank_support.hpp>
#include <sdsl/rank_support_v.hpp>
#include <sdsl/rank_support_v5.hpp>
#include "../src/intersection_high_btrie_log.hpp"
#include <sdsl/vectors.hpp>
#include <sdsl/int_vector.hpp>

#include <limits> // Para obtener los valores máximos y mínimos de los tipos numéricos
#include <numeric> // Para std::accumulate
#include <algorithm> // Para std::sort
#include <cmath>   // Para std::sqrt y std::pow
#include <iomanip> // Para std::setprecision

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

// Define a struct to hold the data you want to write
struct OutputData {
    size_t hsSize;
    double timeRatio;
    size_t intersectionSize;
    size_t clusterCounts[4];
    uint64_t minVal;
    uint64_t maxVal;
    double mean;
    double median;
    double stddev;
    double clusterLimit[3];
};

// Function to calculate the minimum and maximum of a vector of unsigned 64-bit integers
std::pair<uint64_t, uint64_t> calculateMinMax(const std::vector<uint64_t>& vec) {
    if (vec.empty()) {
        return {0, 0};
    }
    uint64_t minVal = vec[0];
    uint64_t maxVal = vec[0];
    for (uint64_t val : vec) {
        if (val < minVal) {
            minVal = val;
        }
        if (val > maxVal) {
            maxVal = val;
        }
    }
    return {minVal, maxVal};
}

// Function to calculate the mean (average) of a vector of unsigned 64-bit integers
double calculateMean(const std::vector<uint64_t>& vec) {
    if (vec.empty()) {
        return 0.0;
    }
    double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    return sum / vec.size();
}

// Function to calculate the median of a vector of unsigned 64-bit integers
double calculateMedian(std::vector<uint64_t> vec) {
    if (vec.empty()) {
        return 0.0;
    }
    std::sort(vec.begin(), vec.end());
    size_t n = vec.size();
    if (n % 2 == 0) {
        // If the size is even, the median is the average of the two central elements
        return (static_cast<double>(vec[n / 2 - 1]) + vec[n / 2]) / 2.0;
    } else {
        // If the size is odd, the median is the central element
        return static_cast<double>(vec[n / 2]);
    }
}

// Function to calculate the standard deviation of a vector of unsigned 64-bit integers
double calculateStandardDeviation(const std::vector<uint64_t>& vec, double mean) {
    if (vec.empty()) {
        return 0.0;
    }
    double squaredSum = 0.0;
    for (uint64_t val : vec) {
        squaredSum += std::pow(val - mean, 2);
    }
    return std::sqrt(squaredSum / vec.size());
}

// Function to count the number of elements in 4 clusters based on the mean
std::vector<size_t> countElementsByMean(const std::vector<uint64_t>& numbers) {
    std::vector<size_t> clusterCounts(4, 0); // Initialize counts to 0
    if (numbers.empty()) {
        return clusterCounts;
    }

    double mean = calculateMean(numbers);
    uint64_t maxElement = *std::max_element(numbers.begin(), numbers.end());
    double limit1 = mean / 2.0;
    double limit2 = mean;
    double limit3 = mean + (maxElement - mean) / 2.0;

    for (uint64_t val : numbers) {
        if (val <= limit1) {
            clusterCounts[0]++;
        } else if (val <= limit2) {
            clusterCounts[1]++;
        } else if (val <= limit3) {
            clusterCounts[2]++;
        } else {
            clusterCounts[3]++;
        }
    }
    return clusterCounts;
}



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
void performIntersectionsHS_High_BTrie_Log( std::ifstream &in_sequences, std::string query_path,
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
        out.open(out_path, std::ios::out|std::ios::binary);
        // out << "elements_per_query,time execution,size_intersection" << std::endl; 
    }

    uint64_t nq = 0;
    long unsigned int total_time = 0;

    for (auto q: queries) {
        vector<HybridStructure> Hs;
        for (auto i: q){
            Hs.push_back(*sequences[i]);            
        }

        vector<uint64_t> intersection_High_BTrie;
        if (Hs.size() <= 16){
            uint64_t time_10 = 0;
            for(int rep = 0; rep < trep; ++rep) {
                auto start = std::chrono::high_resolution_clock::now();
                
                Intersect_HS_High_BTrie_Log<HybridStructure, low_part_size>(Hs, intersection_High_BTrie, intersectBF, runs_encoded, parallel);

                auto end = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
                auto time = elapsed.count();
                total_time += time;
                time_10 += time;
                if (rep != trep-1)
                    intersection_High_BTrie.clear();    
            }
            if (out.is_open()) {
                OutputData dataToWrite;

                dataToWrite.hsSize = Hs.size();
                dataToWrite.timeRatio = static_cast<float>(time_10*1e-3) / trep;
                dataToWrite.intersectionSize = intersection_High_BTrie.size();
                std::vector<size_t> clusterCounts = countElementsByMean(intersection_High_BTrie);
                for (int i = 0; i < 4; ++i) {
                    dataToWrite.clusterCounts[i] = clusterCounts[i];
                }
                std::pair<uint64_t, uint64_t> minMax = calculateMinMax(intersection_High_BTrie);
                dataToWrite.minVal = minMax.first;
                dataToWrite.maxVal = minMax.second;
                dataToWrite.mean = calculateMean(intersection_High_BTrie);
                dataToWrite.median = calculateMedian(intersection_High_BTrie);
                dataToWrite.stddev = calculateStandardDeviation(intersection_High_BTrie, dataToWrite.mean);
                dataToWrite.clusterLimit[0] = dataToWrite.mean / 2.0;
                dataToWrite.clusterLimit[1] = dataToWrite.mean;
                dataToWrite.clusterLimit[2] = dataToWrite.mean + (dataToWrite.maxVal - dataToWrite.mean) / 2.0;

                if (intersection_High_BTrie.size() == 0) {
                    cout << "Min value is 0" << std::endl;
                } else {
                                    // Write the entire struct to the binary file
                    out.write(reinterpret_cast<const char*>(&dataToWrite), sizeof(OutputData));
                }
               
            }
            ++nq;
            if (nq % 1000 == 0) {
                std::cout << nq << " queries processed" << std::endl;
            }
        }
    }

    std::cout << "Number of queries: " << nq << std::endl;
    std::cout <<"Avg time execution: " << (double)(total_time*1e-3)/(nq*trep) << "[ms]" << std::endl;

    out.close();
}
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_scan<1>, 512, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 1>, 1>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_mcl<1>, 512, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_scan<1>, 512, 2>, 2>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);


template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_mcl<1>, 512, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_scan<1>, 512, 4>, 4>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_scan<1>, 512, 8>, 8>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_scan<1>, 512, 12>, 12>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);

template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);
template
void performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_scan<1>, 512, 16>, 16>(std::ifstream &in_sequences, std::string query_path, std::string out_path, bool runs_encoded, bool intersectBF, bool parallel);


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


bool verbose = false;





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
      
        if (std::string(argv[i]) == "--out") {
            ++i;
            output_filename = std::string(argv[i]);
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

    
    std::cout << "outfilename: " << output_filename << endl;



    if(low_p == 1) {

            if (rank == 0 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } else if (rank == 1 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_scan<1>, 512, 1>, 1>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
        
        
    } 
    else if(low_p == 2) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_scan<1>, 512, 2>, 2>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
    }
    else if(low_p == 4) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } else if (rank == 1 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_scan<1>, 512, 4>, 4>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
    } 
    else if(low_p == 8) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } else if (rank == 1 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_scan<1>, 512, 8>, 8>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
        
    } 
    else if(low_p == 12) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_scan<1>, 512, 12>, 12>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
        
    } 
    else if(low_p == 16) {
            if (rank == 0 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_mcl<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 0  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v<1>, select_support_scan<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel); 
            }  
            else if (rank == 2 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_mcl<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }
            else if (rank == 2  && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_flat<rank_support_v5<1>, select_support_scan<1>, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 0 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il<select_support_mcl<1>, 512, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            } 
            else if (rank == 1 && select == 1 ) {
                performIntersectionsHS_High_BTrie_Log<HybridStructure_il< select_support_scan<1>, 512, 16>, 16>(in_sequences, querylog_filename, output_filename, runs, intersectBF, parallel);
            }  
        }
        
    
    
    std::cout << "END" << endl;

    return 0;
}