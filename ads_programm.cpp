#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <cmath>
#include <unordered_map>
#include <algorithm>

using namespace std;

class AbstractRMQ
{
public:
    virtual uint32_t rmq(uint32_t, uint32_t) = 0;
};

class NaiveRMQ : public AbstractRMQ
{
private:
    vector<vector<uint32_t>> rmq_solutions;

public:
    NaiveRMQ(uint32_t n, vector<uint64_t> &v)
    {
        rmq_solutions = vector<vector<uint32_t>>(n - 1);

        for (uint32_t i = 0; i < n - 1; i++)
        {
            uint32_t cur_min_idx = i;
            for (uint32_t j = i + 1; j < n; j++)
            {
                if (v[j] < v[cur_min_idx])
                {
                    cur_min_idx = j;
                }
                rmq_solutions[i].push_back(cur_min_idx);
            }
        }
    }

    uint32_t rmq(uint32_t s, uint32_t e)
    {
        if (s == e)
            return s;
        return rmq_solutions[s][e - s - 1];
    }
};

class LogLinearRMQ : public AbstractRMQ
{
private:
    vector<vector<uint32_t>> rmq_solutions;
    vector<uint64_t> *v;

    bool is_power_of_2(uint32_t n)
    {
        return ((n & (n - 1)) == 0);
    }

public:
    LogLinearRMQ(uint32_t n, vector<uint64_t> &v)
    {
        rmq_solutions = vector<vector<uint32_t>>(n - 1);
        this->v = &v;

        uint32_t idx1, idx2;

        for (uint32_t l = 1; l <= log2(n); l++)
        {
            for (uint32_t s = 0; s + pow(2, l) - 1 < n; s++)
            {
                idx1 = (l > 1) ? rmq_solutions[s][l - 2] : s;
                idx2 = (l > 1) ? rmq_solutions[s + pow(2, l - 1)][l - 2] : s + pow(2, l - 1);
                if (v[idx1] < v[idx2])
                {
                    rmq_solutions[s].push_back(idx1);
                }
                else
                {
                    rmq_solutions[s].push_back(idx2);
                }
            }
        }
    }

    uint32_t rmq(uint32_t s, uint32_t e)
    {
        uint32_t l = log2(e - s + 1);

        if (is_power_of_2(e - s + 1))
        {
            return (l > 0) ? rmq_solutions[s][l - 1] : s;
        }

        uint32_t idx1 = (l > 0) ? rmq_solutions[s][l - 1] : s;
        uint32_t idx2 = (l > 0) ? rmq_solutions[e - pow(2, l) + 1][l - 1] : e - pow(2, l) + 1;

        if ((*v)[idx1] < (*v)[idx2])
        {
            return idx1;
        }
        else
        {
            return idx2;
        }
    }
};

struct CartesianTree
{
    uint32_t min_idx;
    CartesianTree *parent = nullptr;
    CartesianTree *left_child = nullptr;
    CartesianTree *right_child = nullptr;
    vector<bool> representation;

    bool is_root() { return parent == nullptr; }
    // vector<bool> get_unbalanced_binary_dfuds_representation()
    // {
    //     vector<bool> dfuds;
    //     if (left_child != nullptr)
    //     {
    //         dfuds.push_back(0);
    //     }
    //     if (right_child != nullptr)
    //     {
    //         dfuds.push_back(0);
    //     }
    //     dfuds.push_back(1);
    //     if (left_child != nullptr)
    //     {
    //         vector<bool> left_dfuds = left_child->get_unbalanced_binary_dfuds_representation();
    //         dfuds.insert(dfuds.end(), left_dfuds.begin(), left_dfuds.end());
    //     }
    //     if (right_child != nullptr)
    //     {
    //         vector<bool> right_dfuds = right_child->get_unbalanced_binary_dfuds_representation();
    //         dfuds.insert(dfuds.end(), right_dfuds.begin(), right_dfuds.end());
    //     }

    //     return dfuds;
    // }
    // uint32_t rmq(uint32_t s, uint32_t e) {
    //     if ((s <= min_idx) && (min_idx <= e)) {
    //         return min_idx;
    //     } else if (s > min_idx) {
    //         return left_child->rmq(s, e);
    //     } else if (e < min_idx) {
    //         return right_child->rmq(s, e);
    //     } else {
    //         throw runtime_error("Invalid parameters for rmq:" + to_string(s) + "-" + to_string(e));
    //     }
    // }
};

template <typename T>
void printVector(const std::vector<T>& vec) {
    for(T i : vec) {
        std::cout << i << " ";
    }
    std::cout << "\n";
}

CartesianTree construct_cartesian_tree(vector<uint64_t> &v, uint32_t start_idx, uint32_t end_idx, uint32_t block_size)
{      
    vector<bool> repr;
    CartesianTree* c_tree = new CartesianTree();
    c_tree->min_idx = 0;

    CartesianTree* cur_c_tree = c_tree;
    for (uint32_t i = start_idx + 1; i < end_idx; i++)
    {   
        while ((!cur_c_tree->is_root()) && (v[i] < v[start_idx + cur_c_tree->min_idx]))
        {
            repr.push_back(1);
            cur_c_tree = cur_c_tree->parent;
        }

        CartesianTree *new_c_tree = new CartesianTree();
        new_c_tree->min_idx = i - start_idx;

        if (cur_c_tree->is_root() && (v[i] < v[start_idx + cur_c_tree->min_idx]))
        {
            repr.push_back(1);
            new_c_tree->left_child = cur_c_tree;
            cur_c_tree->parent = new_c_tree;
            c_tree = new_c_tree;
        }
        else
        {
            new_c_tree->parent = cur_c_tree;
            new_c_tree->left_child = cur_c_tree->right_child;
            if (new_c_tree->left_child != nullptr) {
                new_c_tree->left_child->parent = new_c_tree;
            }
            cur_c_tree->right_child = new_c_tree;
        }

        repr.push_back(0);

        cur_c_tree = new_c_tree;
    }

    for (uint32_t i = end_idx - start_idx; i < block_size; i++) {
        repr.push_back(0);
    }

    c_tree->representation = repr;
    return *c_tree;
}

class LinearRMQ : public AbstractRMQ
{
private:
    uint32_t block_size;
    vector<uint64_t> *v;

    vector<uint64_t> min_within_block;
    vector<uint32_t> min_idx_within_block;
    LogLinearRMQ *query_spanning_block_rmq_ds;

    vector<vector<bool>> c_trees_dfuds;
    unordered_map<vector<bool>, vector<vector<uint32_t>>> c_tree_start_end_rmqs;

    void construct_c_tree_start_end_rmqs(uint32_t n) {
        vector<uint64_t> v(n);
    
        // Initialize vector with values from 0 to n-1
        for(uint32_t i = 0; i < n; i++) {
            v[i] = i;
        }

        do {
            LogLinearRMQ rmq_ds = LogLinearRMQ(n, v);
            CartesianTree c_tree = construct_cartesian_tree(v, 0, n, n);
            // vector<bool> dfuds_repr = c_tree.get_unbalanced_binary_dfuds_representation();
            vector<bool> dfuds_repr = c_tree.representation;

            c_tree_start_end_rmqs[dfuds_repr] = vector<vector<uint32_t>>(n, vector<uint32_t>(n));

            for (uint32_t i = 0; i < n; i++) {
                for (uint32_t j = i; j < n; j++) {
                    c_tree_start_end_rmqs[dfuds_repr][i][j] = rmq_ds.rmq(i, j);
                }
            }
        } while (next_permutation(v.begin(), v.end()));
    }

public:
    LinearRMQ(uint32_t n, vector<uint64_t> &v)
    {
        this->v = &v;
        block_size = log2(n) / 4;

        construct_c_tree_start_end_rmqs(block_size);

        for (uint32_t i = 0; i < ceil((double) n / block_size); i++)
        {
            uint32_t start_idx = block_size * i;
            uint32_t end_idx = min(block_size * (i + 1), n); 

            uint64_t min_this_block = v[start_idx];
            uint32_t min_idx_this_block = start_idx;

            for (uint32_t j = start_idx + 1; j < end_idx; j++)
            {
                if (v[j] < min_this_block)
                {
                    min_this_block = v[j];
                    min_idx_this_block = j;
                }
            }

            min_within_block.push_back(min_this_block);
            min_idx_within_block.push_back(min_idx_this_block);

            CartesianTree block_c_tree = construct_cartesian_tree(v, start_idx, end_idx, block_size);
            c_trees_dfuds.push_back(block_c_tree.representation);
        }

        query_spanning_block_rmq_ds = new LogLinearRMQ(n, min_within_block);
    }

    uint32_t spanning_block_rmq(uint32_t s_block, uint32_t e_block)
    {
        uint32_t min_block_idx = query_spanning_block_rmq_ds->rmq(s_block, e_block);
        return min_idx_within_block[min_block_idx];
    }

    uint32_t within_block_rmq(uint32_t block_idx, uint32_t s, uint32_t e) {
        return (block_idx * block_size) + c_tree_start_end_rmqs[c_trees_dfuds[block_idx]][s][e];
    }

    uint32_t rmq(uint32_t s, uint32_t e)
    {
        uint32_t s_block = (double) s / block_size;
        uint32_t e_block = (double) e / block_size;

        if (s_block == e_block) {
            return within_block_rmq(s_block, s % block_size, e % block_size);
        } else if ((s % block_size == 0) && ((e % block_size == block_size - 1) || (e == v->size() - 1))) {
            return spanning_block_rmq(s_block, e_block);
        } else if (s % block_size == 0) {
            uint32_t within_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block, e_block - 1);
            return ((*v)[within_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_block_min_idx : spanning_block_min_idx;
        } else if ((e % block_size == block_size - 1) || (e == v->size() - 1)) {
            uint32_t within_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block + 1, e_block);
            return ((*v)[within_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_block_min_idx : spanning_block_min_idx;
        } else if (s_block + 1 == e_block) {
            uint32_t within_s_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t within_e_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            return ((*v)[within_s_block_min_idx] < (*v)[within_e_block_min_idx]) ? within_s_block_min_idx : within_e_block_min_idx;
        } else {
            uint32_t within_s_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t within_e_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block + 1, e_block - 1);

            if ((*v)[within_s_block_min_idx] < (*v)[within_e_block_min_idx]) {
                return ((*v)[within_s_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_s_block_min_idx : spanning_block_min_idx;
            } else {
                return ((*v)[within_e_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_e_block_min_idx : spanning_block_min_idx;
            }

        }
    }
};

void run_naive_rmq(ifstream &input_file, ofstream &output_file)
{
    vector<uint64_t> v;
    uint32_t n;

    input_file >> n;

    uint64_t elem;
    for (uint32_t i = 0; i < n; i++)
    {
        input_file >> elem;
        v.push_back(elem);
    }

    NaiveRMQ rmq_ds = NaiveRMQ(n, v);
    v.clear();

    uint32_t s, e;
    string line;

    while (input_file >> line)
    {
        size_t comma_idx = line.find(',');
        uint32_t s = stoi(line.substr(0, comma_idx));
        uint32_t e = stoi(line.substr(comma_idx + 1));

        output_file << rmq_ds.rmq(s, e) << endl;
    }
}

void run_loglinear_rmq(ifstream &input_file, ofstream &output_file)
{
    vector<uint64_t> v;
    uint32_t n;

    input_file >> n;

    uint64_t elem;
    for (uint32_t i = 0; i < n; i++)
    {
        input_file >> elem;
        v.push_back(elem);
    }

    LogLinearRMQ rmq_ds = LogLinearRMQ(n, v);

    uint32_t s, e;
    string line;

    while (input_file >> line)
    {
        size_t comma_idx = line.find(',');
        uint32_t s = stoi(line.substr(0, comma_idx));
        uint32_t e = stoi(line.substr(comma_idx + 1));

        output_file << rmq_ds.rmq(s, e) << endl;
    }
}

void run_linear_rmq(ifstream &input_file, ofstream &output_file)
{
    vector<uint64_t> v;
    uint32_t n;

    input_file >> n;

    uint64_t elem;
    for (uint32_t i = 0; i < n; i++)
    {
        input_file >> elem;
        v.push_back(elem);
    }

    AbstractRMQ* rmq_ds;

    if (((uint32_t) log2(n) / 4) == 0) {
        rmq_ds = new LogLinearRMQ(n, v);
    } else {
        rmq_ds = new LinearRMQ(n, v);
    }


    uint32_t s, e;
    string line;

    while (input_file >> line)
    {
        size_t comma_idx = line.find(',');
        uint32_t s = stoi(line.substr(0, comma_idx));
        uint32_t e = stoi(line.substr(comma_idx + 1));

        output_file << rmq_ds->rmq(s, e) << endl;
    }
}


void test()
{
    vector<uint64_t> v{0, 2, 1, 3};
    CartesianTree c_tree = construct_cartesian_tree(v, 0, 4, 4);
    printVector(c_tree.representation);
}

int main(int argc, char *argv[])
{
    // test();
    // exit(0);
    ifstream input_file(argv[2]);
    if (!input_file.is_open())
    {
        throw runtime_error("Unable to open input file "s + argv[2]);
    }
    ofstream output_file(argv[3]);
    if (!output_file.is_open())
    {
        throw runtime_error("Unable to open output file "s + argv[3]);
    }

    if (argv[1] == "pd"s)
    {
    }
    else if (argv[1] == "rmq"s)
    {
        run_linear_rmq(input_file, output_file);
    }
    else
    {
        throw invalid_argument("The first parameter must either be 'pd' or 'rmq'.");
    }
}

// int main() {
//     int n;
//     std::cout << "Enter a number: ";
//     std::cin >> n;
    
//     std::vector<int> vec(n);
    
//     // Initialize vector with values from 0 to n-1
//     for(int i = 0; i < n; i++) {
//         vec[i] = i;
//     }
    
//     // Print the initial vector
//     printVector(vec);

//     // Generate all possible permutations
//     while(std::next_permutation(vec.begin(), vec.end())) {
//         printVector(vec);
//     }
    
//     return 0;
// }