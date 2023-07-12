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
    virtual int rmq(uint32_t ,uint32_t) = 0;
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

        int rmq(uint32_t s, uint32_t e) {
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
            rmq_solutions = vector<vector<uint32_t>>(n);
            this->v = &v;

            for (uint32_t s = 0; s < n; s++) {
                rmq_solutions[s].push_back(s);
            }
            
            uint32_t idx1, idx2;

            for (uint32_t l = 1; l <= log2(n); l++) {
                for (uint32_t s = 0; s + pow(2, l) - 1 < n; s++) {
                    idx1 = rmq_solutions[s][l - 1];
                    idx2 = rmq_solutions[s + pow(2, l - 1)][l - 1];
                    if (v[idx1] < v[idx2]) {
                        rmq_solutions[s].push_back(idx1);
                    } else {
                        rmq_solutions[s].push_back(idx2);
                    }
                }
            }
        }

        int rmq(uint32_t s, uint32_t e) {
            if (is_power_of_2(e - s + 1)) {
                uint32_t l = log2(e - s + 1);
                return rmq_solutions[s][l];
            }
            uint32_t l = floor(log2(e - s - 1));

            uint32_t idx1 = rmq_solutions[s][l];
            uint32_t idx2 = rmq_solutions[e - pow(2, l) + 1][l];

            if ((*v)[idx1] < (*v)[idx2]) {
                return idx1;
            } else {
                return idx2;
            }
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