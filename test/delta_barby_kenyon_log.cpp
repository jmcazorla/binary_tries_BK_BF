#include <iostream>
#include <vector>
#include <fstream>
#include "../src/barbay_and_kenyon.hpp"
#include <sdsl/vectors.hpp>


using namespace std;
using namespace sdsl;



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
        if(query.size() > 1 ) 
            queries->push_back(query);
    }
    in.close();
    return queries;
}

std::vector<uint64_t>* read_inv_list(std::ifstream &input_stream, uint32_t n) {

    std::vector <uint64_t>* il = new std::vector<uint64_t>();
    il -> reserve(n);
    for (uint32_t i = 0; i < n; ++i) {
        uint32_t x;
        input_stream.read(reinterpret_cast<char *>(&x), 4);
        il -> push_back((uint64_t)x);
    }
    return il;
}


void performDeltaBarbayKenyon( std::string sequences_path, std::string query_path,
                           std::string out_path, uint64_t min_size, uint64_t max_size) {
    
    
    std::ifstream input_stream(sequences_path, std::ios::binary | std::ios::in);

    if (!input_stream.is_open() ) {
        cout << "Cant't open file: " << sequences_path << endl;
    }

    uint32_t _1, u;
    // Read universe of collection
    input_stream.read(reinterpret_cast<char *>(&_1), sizeof(_1));
    input_stream.read(reinterpret_cast<char *>(&u), sizeof(u));

    uint32_t nSets = 0;
    while (true) {
        uint32_t n;
        input_stream.read(reinterpret_cast<char *>(&n), 4);
        if (input_stream.eof()){
            break;
        }
        if (n > min_size && (n < max_size || max_size == 0)) {
            nSets++;
        }
        input_stream.seekg(4*n, std::ios::cur);
    }
    // Set file pointer in first set
    input_stream.clear();
    input_stream.seekg(2*sizeof(u), std::ios::beg);
    // Write universe in out file
    std::ofstream out;
    if (out_path != "") {
        out.open(out_path, std::ios::out);
        out << "universo, number_of_sets, min_size, max_size," << std::endl; 
        out << u << ","<< nSets << "," << min_size<< "," << max_size << std::endl; 
        out << "query, sets_per_query, elements_per_query, delta" << std::endl; 
        out << "elements_per_sets 1, elements_persets 2, ...." << std::endl; 
        out << "number_of_set_in_query 1, number_of_set_in_query 2, ...." << std::endl; 
    }


    std::cout << "Universo: " << u << std::endl;
    std::cout << "Number of sets: " << nSets << std::endl;


    vector<vector<uint32_t>>* queries;
    queries   = loadQueryLog(query_path, nSets);
    cout << "Queries loaded succefully: "  << queries->size() << std::endl<< std::endl;

    vector<vector<uint64_t>> *sequences;
    uint64_t n_il = 0;
    uint64_t nq = 0;

    for (int i = 0; i < queries -> size(); ++i) {
        sequences = new vector<vector<uint64_t>>();
        vector<uint32_t> query = (*queries)[i];
        std::sort(query.begin(), query.end());
        uint32_t count_set_query = 0;
        uint64_t total_set_elements = 0;
        uint32_t read_set_size = -1;

        while (true && query.size() > count_set_query) {
            uint32_t n;
            if (!input_stream.read(reinterpret_cast<char*>(&n), sizeof(uint32_t))) {
                break; 
            }            
            read_set_size ++;
           if (n > min_size && (n < max_size || max_size == 0)) {
                    
                if ( query[count_set_query] == read_set_size){ 
                    total_set_elements += n;
                    count_set_query++;

                    std::vector <uint64_t> *il = read_inv_list(input_stream, n);
                    sequences->push_back(*il);
                    n_il++;
                    delete il;
                } else {
                    input_stream.seekg(4*n, std::ios::cur);
                }
            } else {
                input_stream.seekg(4*n, std::ios::cur);
            }        
        }

        if (count_set_query <= 16 && count_set_query == query.size()){
            uint64_t delta = 0;
            deltaBarbayKenyon((*sequences), count_set_query, delta);
            if (out.is_open()) {
                out << i << ","<< count_set_query << ","<< total_set_elements << "," << delta << std::endl;
                for (int k = 0; k < sequences -> size(); ++k) {
                    k+1 == sequences -> size() ? out<<(*sequences)[k].size()<< std::endl : out<<(*sequences)[k].size()<< ",";
                }
                for (int k = 0; k < query.size(); ++k) {
                    k+1 == query.size() ? out<<query[k]<< std::endl : out<<query[k]<< ",";
                }
                out << std::endl;	
            }
            ++nq;
            if (nq % 1000 == 0) {
                std::cout << nq << " queries processed" << std::endl;
            }

        } else {
             if (out.is_open()) {
                out << std::endl;
             }
        }
       

        delete sequences;

        // Set file pointer in first set
        input_stream.clear();
        input_stream.seekg(2*sizeof(u), std::ios::beg);
        // Write universe in out file
    }
    std::cout << "Total inverted lists: " << n_il << std::endl;



    input_stream.close();
    out.close();

    delete queries;
}

void performQueriesSets( std::string sequences_path, std::string query_path,
                           std::string out_path, uint64_t min_size, uint64_t max_size) {
    
    
    std::ifstream input_stream(sequences_path, std::ios::binary | std::ios::in);

    if (!input_stream.is_open() ) {
        cout << "Cant't open file: " << sequences_path << endl;
    }

    uint32_t _1, u;
    // Read universe of collection
    input_stream.read(reinterpret_cast<char *>(&_1), sizeof(_1));
    input_stream.read(reinterpret_cast<char *>(&u), sizeof(u));

    uint32_t nSets = 0;
    while (true) {
        uint32_t n;
        input_stream.read(reinterpret_cast<char *>(&n), 4);
        if (input_stream.eof()){
            break;
        }
        if (n > min_size && (n < max_size || max_size == 0)) {
            nSets++;
        }
        input_stream.seekg(4*n, std::ios::cur);
    }
    // Set file pointer in first set
    input_stream.clear();
    input_stream.seekg(2*sizeof(u), std::ios::beg);
    // Write universe in out file
    std::ofstream out;
    if (out_path != "") {
        out.open(out_path, std::ios::out);
        out << "universo, number_of_sets, min_size, max_size," << std::endl; 
        out << u << ","<< nSets << "," << min_size<< "," << max_size << std::endl; 
        out << "query, sets_per_query, elements_per_query" << std::endl; 
        out << "sets 1, sets 2, ...." << std::endl; 
        out << "number_of_set_in_query 1, number_of_set_in_query 2, ...." << std::endl; 
    }


    std::cout << "Universo: " << u << std::endl;
    std::cout << "Number of sets: " << nSets << std::endl;


    vector<vector<uint32_t>>* queries;
    queries   = loadQueryLog(query_path, nSets);
    cout << "Queries loaded succefully: "  << queries->size() << std::endl<< std::endl;

    vector<vector<uint64_t>> *sequences;
    uint64_t n_il = 0;
    uint64_t nq = 0;

    for (int i = 0; i < queries->size(); ++i) {
        sequences = new vector<vector<uint64_t>>();
        vector<uint32_t> query = (*queries)[i];
        std::sort(query.begin(), query.end());
        uint32_t count_set_query = 0;
        uint64_t total_set_elements = 0;
        uint32_t read_set_size = -1;

        while (true && query.size() > count_set_query) {
            uint32_t n;
            if (!input_stream.read(reinterpret_cast<char*>(&n), sizeof(uint32_t))) {
                break; 
            }
            

           if (n > min_size && (n < max_size || max_size == 0)) {
                read_set_size ++;
                    
                if ( query[count_set_query] == read_set_size){ 
                    total_set_elements += n;
                    count_set_query++;
                    std::vector <uint64_t> *il = read_inv_list(input_stream, n);
                    sequences->push_back(*il);
                    n_il++;
                    delete il;
                } else {
                    input_stream.seekg(4*n, std::ios::cur);
                }

            } else {
                input_stream.seekg(4*n, std::ios::cur);
            }
        
        }
        if (count_set_query <= 16 && count_set_query == query.size()){
            uint64_t delta = 0;
            if (out.is_open()) {
                out << i << ","<< count_set_query << ","<< total_set_elements  << std::endl;
                for (int k = 0; k < sequences->size(); ++k) {
                    k+1 == sequences -> size() ? out<<(*sequences)[k].size()<< std::endl : out<<(*sequences)[k].size()<< ",";
                }
                for (int k = 0; k < query.size(); ++k) {
                    k+1 == query.size() ? out<<query[k]<< std::endl : out<<query[k]<< ",";
                }
            }
            ++nq;
            if (nq % 1000 == 0) {
                std::cout << nq << " queries processed" << std::endl;
            }
        }

        delete sequences;

        // Set file pointer in first set
        input_stream.clear();
        input_stream.seekg(2*sizeof(u), std::ios::beg);
        // Write universe in out file
    }
    std::cout << "Total inverted lists: " << n_il << std::endl;



    input_stream.close();
    out.close();

    delete queries;
}


int main(int argc, char const *argv[]) {
    int mandatory = 5;
    // (*) mandatory parameters
    if (argc < mandatory){
        std::cout   << "collection filename " // (*)
                        "query log "          // (*)
                        "[--min_size m] "     // (*)
                        "[--max_size m] "     // (*)
                        "[--is_delta m] " 
                        "[--out results_filename]"
                    <<
        std::endl;
        return 1;
    }

  
    uint64_t min_size = 0;
    uint64_t max_size = 0;
    uint8_t is_delta = 1;
    std::string sequences_filename = std::string(argv[1]);
    std::string querylog_filename   = std::string(argv[2]);
    std::string output_filename = "";
    
    
    for (int i = 1; i < argc; ++i){
        if (std::string(argv[i]) == "--min_size") {
            ++i;
            min_size = std::stoull(argv[i]);
        }
        if (std::string(argv[i]) == "--max_size") {
            ++i;
            max_size = std::stoull(argv[i]);
        }
        if (std::string(argv[i]) == "--is_delta") {
            ++i;
            is_delta = std::stoull(argv[i]);
        }
        if (std::string(argv[i]) == "--out") {
            ++i;
            output_filename = std::string(argv[i]);
        }
    }
    
    std::cout << "Min size: " << min_size << std::endl;
    std::cout << "Max size: " << max_size << std::endl;
    std::cout << "outfilename: " << output_filename << endl;


    is_delta == 1 ? performDeltaBarbayKenyon(sequences_filename,querylog_filename,output_filename,min_size, max_size)
                    : performQueriesSets(sequences_filename,querylog_filename,output_filename,min_size, max_size);



    return 0;
}