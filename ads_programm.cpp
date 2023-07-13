#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <cmath>

using namespace std;

class AbstractRMQ
{
    public:
        virtual uint32_t rmq(uint32_t ,uint32_t) = 0;
};

class NaiveRMQ : public AbstractRMQ {
    private:
        vector<vector<uint32_t>> rmq_solutions;
    public:
        NaiveRMQ(uint32_t n, vector<uint64_t> &v) {
            rmq_solutions = vector<vector<uint32_t>>(n - 1);

            for (uint32_t i = 0; i < n - 1; i++) {
                uint32_t cur_min_idx = i;
                for (uint32_t j = i + 1; j < n; j++) {
                    if (v[j] < v[cur_min_idx]) {
                        cur_min_idx = j;
                    }
                    rmq_solutions[i].push_back(cur_min_idx);
                }
            }

        }

        uint32_t rmq(uint32_t s, uint32_t e) {
            if (s == e) return s;
            return rmq_solutions[s][e - s - 1];
        }
};

class LogLinearRMQ : public AbstractRMQ {
    private:
        vector<vector<uint32_t>> rmq_solutions;
        vector<uint64_t> *v;

        bool is_power_of_2(uint32_t n) {
            return ((n & (n - 1)) == 0);
        }
    public:
        LogLinearRMQ(uint32_t n, vector<uint64_t> &v) {
            rmq_solutions = vector<vector<uint32_t>>(n - 1);
            this->v = &v;
            
            uint32_t idx1, idx2;

            for (uint32_t l = 1; l <= log2(n); l++) {
                for (uint32_t s = 0; s + pow(2, l) - 1 < n; s++) {
                    idx1 = (l > 1) ? rmq_solutions[s][l - 2] : s;
                    idx2 = (l > 1) ? rmq_solutions[s + pow(2, l - 1)][l - 2] : s + pow(2, l - 1);
                    if (v[idx1] < v[idx2]) {
                        rmq_solutions[s].push_back(idx1);
                    } else {
                        rmq_solutions[s].push_back(idx2);
                    }
                }
            }

        }

        uint32_t rmq(uint32_t s, uint32_t e) {
            uint32_t l = log2(e - s + 1);

            if (is_power_of_2(e - s + 1)) {
                return (l > 0) ? rmq_solutions[s][l - 1] : s;
            }

            uint32_t idx1 = (l > 0) ? rmq_solutions[s][l - 1] : s;
            uint32_t idx2 = (l > 0) ? rmq_solutions[e - pow(2, l) + 1][l - 1]: e - pow(2, l) + 1;

            if ((*v)[idx1] < (*v)[idx2]) {
                return idx1;
            } else {
                return idx2;
            }
        }
};

class LinearRMQ : public AbstractRMQ {
    private:
        vector<uint64_t> min_within_block;
        vector<uint32_t> min_idx_within_block;
        LogLinearRMQ* query_spanning_block_rmq_ds;
    public:
        LinearRMQ(uint32_t n, vector<uint64_t> &v) {
            double s = log(n) / 4;

            for (uint32_t i = 0; i < ceil(n / s); i++) {
                uint64_t min_this_block = v[0];
                uint32_t min_idx_this_block = 0;
                for (uint32_t j = 1; j < s; j++) {
                    if (v[j] < min_this_block) {
                        min_this_block = v[j];
                        min_idx_this_block = j;
                    }
                }

                min_within_block.push_back(min_this_block);
                min_idx_within_block.push_back(min_idx_this_block);
            }

            query_spanning_block_rmq_ds = new LogLinearRMQ(n, v);
        }

        uint32_t spanning_block_rmq(uint32_t block_idx1, uint32_t block_idx2) {
            uint32_t min_block_idx = query_spanning_block_rmq_ds->rmq(block_idx1, block_idx2);
            return min_idx_within_block[min_block_idx];
        }

        uint32_t rmq(uint32_t s, uint32_t e) {
            // TODO
            return 0;
        }
};


void run_naive_rmq(ifstream &input_file, ofstream &output_file) {
    vector<uint64_t> v;
    uint32_t n;

    input_file >> n;

    uint64_t elem;
    for (uint32_t i = 0; i < n; i++) {
        input_file >> elem;
        v.push_back(elem);
    }

    NaiveRMQ rmq_ds = NaiveRMQ(n, v);
    v.clear();

    uint32_t s, e;
    string line;

    while (input_file >> line) {
        size_t comma_idx = line.find(','); 
        uint32_t s = stoi(line.substr(0, comma_idx));
        uint32_t e = stoi(line.substr(comma_idx + 1));

        output_file << rmq_ds.rmq(s, e) << endl;
    }
}

void run_loglinear_rmq(ifstream &input_file, ofstream &output_file) {
    vector<uint64_t> v;
    uint32_t n;

    input_file >> n;

    uint64_t elem;
    for (uint32_t i = 0; i < n; i++) {
        input_file >> elem;
        v.push_back(elem);
    }

    LogLinearRMQ rmq_ds = LogLinearRMQ(n, v);

    uint32_t s, e;
    string line;

    while (input_file >> line) {
        size_t comma_idx = line.find(','); 
        uint32_t s = stoi(line.substr(0, comma_idx));
        uint32_t e = stoi(line.substr(comma_idx + 1));

        output_file << rmq_ds.rmq(s, e) << endl;
    }
}

int main(int argc, char* argv[]) {

    ifstream input_file (argv[2]);
    if (!input_file.is_open()) {
        throw runtime_error("Unable to open input file "s + argv[2]);
    }
    ofstream output_file (argv[3]);
    if (!output_file.is_open()) {
        throw runtime_error("Unable to open output file "s + argv[3]);
    }

    if (argv[1] == "pd"s) {

    } else if (argv[1] == "rmq"s) {
        run_loglinear_rmq(input_file, output_file);
    } else {
        throw invalid_argument("The first parameter must either be 'pd' or 'rmq'.");
    }
}