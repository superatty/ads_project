#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <stdexcept>

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

            for (uint32_t i = 0; i < v.size() - 1; i++) {
                uint32_t cur_min_idx = i;
                for (uint32_t j = i + 1; j < v.size(); j++) {
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

    } else {
        throw invalid_argument("The first parameter must either be 'pd' or 'rmq'.");
    }
}