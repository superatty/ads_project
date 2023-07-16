#include <vector>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include <tuple>

using namespace std;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

/**
 * Abtract class for
 *
 */
class AbstractRMQ
{
public:
    virtual uint32_t rmq(uint32_t, uint32_t) = 0;
    virtual uint64_t size_in_bits() = 0;
};

class NaiveRMQ : public AbstractRMQ
{
private:
    vector<vector<uint32_t>> rmq_solutions;

public:
    NaiveRMQ(vector<uint64_t> &v)
    {
        uint32_t n = v.size();
        rmq_solutions = vector<vector<uint32_t>>(n - 1);

        for (uint32_t i = 0; i < n - 1; i++)
        {
            uint32_t cur_min_idx = i;
            for (uint32_t j = i + 1; j < n; j++)
            {
                if (v[j] < v[cur_min_idx])
                    cur_min_idx = j;
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

    uint64_t size_in_bits()
    {
        uint64_t size_in_bits = 0;
        for (auto &rmq_sol : rmq_solutions)
            size_in_bits += rmq_sol.size() * 32;

        return size_in_bits;
    }
};

class LogLinearRMQ : public AbstractRMQ
{
private:
    vector<vector<uint32_t>> rmq_solutions;
    vector<uint64_t> *v;

public:
    LogLinearRMQ(vector<uint64_t> &v)
    {
        uint32_t n = v.size();
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
                    rmq_solutions[s].push_back(idx1);
                else
                    rmq_solutions[s].push_back(idx2);
            }
        }
    }

    uint32_t rmq(uint32_t s, uint32_t e)
    {
        if (s == e)
            return s;
        uint32_t l = log2(e - s);

        uint32_t idx1 = (l > 0) ? rmq_solutions[s][l - 1] : s;
        uint32_t idx2 = (l > 0) ? rmq_solutions[e - pow(2, l) + 1][l - 1] : e - pow(2, l) + 1;

        if ((*v)[idx1] < (*v)[idx2])
            return idx1;
        else
            return idx2;
    }

    uint64_t size_in_bits()
    {
        uint64_t size_in_bits = (*v).size() * 64;

        for (auto &rmq_sol : rmq_solutions)
            size_in_bits += rmq_sol.size() * 32;

        return size_in_bits;
    }
};

vector<bool> construct_c_tree(vector<uint64_t> &v, uint32_t start_idx, uint32_t end_idx, uint32_t block_size)
{
    vector<bool> repr; // TODO would allocating space be better or worse?
    vector<uint64_t> stack;

    stack.push_back(v[start_idx]);

    for (size_t i = start_idx + 1; i < end_idx; i++)
    {
        while ((!stack.empty()) && (v[i] < stack.back()))
        {
            repr.push_back(1);
            stack.pop_back();
        }

        repr.push_back(0);
        stack.push_back(v[i]);
    }

    for (uint32_t i = end_idx - start_idx; i < block_size; i++)
        repr.push_back(0);

    return repr;
}

class LinearRMQ : public AbstractRMQ
{
private:
    uint32_t block_size;
    vector<uint64_t> *v;

    vector<uint64_t> min_within_block;
    vector<uint32_t> min_idx_within_block;
    AbstractRMQ *query_spanning_block_rmq_ds;

    vector<vector<bool>> c_trees;
    unordered_map<vector<bool>, vector<vector<uint32_t>>> c_tree_start_end_rmqs;

    void construct_c_tree_start_end_rmqs(uint32_t n)
    {
        vector<uint64_t> v(n);

        // Initialize vector with values from 0 to n-1
        for (uint32_t i = 0; i < n; i++)
            v[i] = i;

        do
        {
            AbstractRMQ *rmq_ds = new NaiveRMQ(v);
            vector<bool> c_tree = construct_c_tree(v, 0, n, n);

            if (c_tree_start_end_rmqs.find(c_tree) != c_tree_start_end_rmqs.end())
                continue;

            c_tree_start_end_rmqs[c_tree] = vector<vector<uint32_t>>(n);

            for (uint32_t i = 0; i < n - 1; i++)
            {
                c_tree_start_end_rmqs[c_tree][i] = vector<uint32_t>(n - i);
                for (uint32_t j = i + 1; j < n; j++)
                    c_tree_start_end_rmqs[c_tree][i][j - 1] = rmq_ds->rmq(i, j);
            }
        } while (next_permutation(v.begin(), v.end()));
    }

public:
    LinearRMQ(vector<uint64_t> &v)
    {
        this->v = &v;
        uint32_t n = v.size();
        block_size = ceil(log2(n) / 4); // taking the ceiling leads to faster execution with reduced space
        uint32_t number_of_blocks = ceil((double)n / block_size);

        min_within_block = vector<uint64_t>(number_of_blocks);
        min_idx_within_block = vector<uint32_t>(number_of_blocks);
        c_trees = vector<vector<bool>>(number_of_blocks);

        construct_c_tree_start_end_rmqs(block_size);

        for (uint32_t i = 0; i < number_of_blocks; i++)
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

            min_within_block[i] = min_this_block;
            min_idx_within_block[i] = min_idx_this_block;

            vector<bool> block_c_tree = construct_c_tree(v, start_idx, end_idx, block_size);
            c_trees[i] = block_c_tree;
        }

        if (number_of_blocks > 100)
            query_spanning_block_rmq_ds = new LinearRMQ(min_within_block);
        else
            query_spanning_block_rmq_ds = new LogLinearRMQ(min_within_block);
    }

    uint32_t spanning_block_rmq(uint32_t s_block, uint32_t e_block)
    {
        uint32_t min_block_idx = query_spanning_block_rmq_ds->rmq(s_block, e_block);
        return min_idx_within_block[min_block_idx];
    }

    uint32_t within_block_rmq(uint32_t block_idx, uint32_t s, uint32_t e)
    {
        if (s == e)
            return (block_idx * block_size) + s;
        return (block_idx * block_size) + c_tree_start_end_rmqs[c_trees[block_idx]][s][e - 1];
    }

    uint32_t rmq(uint32_t s, uint32_t e)
    {
        uint32_t s_block = (double)s / block_size;
        uint32_t e_block = (double)e / block_size;

        if (s_block == e_block)
        {
            return within_block_rmq(s_block, s % block_size, e % block_size);
        }
        else if ((s % block_size == 0) && ((e % block_size == block_size - 1) || (e == v->size() - 1)))
        {
            return spanning_block_rmq(s_block, e_block);
        }
        else if (s % block_size == 0)
        {
            uint32_t within_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block, e_block - 1);
            return ((*v)[within_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_block_min_idx : spanning_block_min_idx;
        }
        else if ((e % block_size == block_size - 1) || (e == v->size() - 1))
        {
            uint32_t within_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block + 1, e_block);
            return ((*v)[within_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_block_min_idx : spanning_block_min_idx;
        }
        else if (s_block + 1 == e_block)
        {
            uint32_t within_s_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t within_e_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            return ((*v)[within_s_block_min_idx] < (*v)[within_e_block_min_idx]) ? within_s_block_min_idx : within_e_block_min_idx;
        }
        else
        {
            uint32_t within_s_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t within_e_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block + 1, e_block - 1);

            if ((*v)[within_s_block_min_idx] < (*v)[within_e_block_min_idx])
            {
                return ((*v)[within_s_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_s_block_min_idx : spanning_block_min_idx;
            }
            else
            {
                return ((*v)[within_e_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_e_block_min_idx : spanning_block_min_idx;
            }
        }
    }

    uint64_t size_in_bits()
    {
        uint64_t size_in_bits = 32;
        size_in_bits += (*v).size() * 64;
        size_in_bits += min_within_block.size() * 64;
        size_in_bits += min_idx_within_block.size() * 32;
        size_in_bits += query_spanning_block_rmq_ds->size_in_bits();

        for (auto c_tree : c_trees)
            size_in_bits += c_tree.size();

        for (auto &start_end_rmq : c_tree_start_end_rmqs)
            for (auto &end_rmq : start_end_rmq.second)
                size_in_bits += end_rmq.size() * 32;

        return size_in_bits;
    }
};

class AbstractBV
{
public:
    virtual uint32_t select0(uint32_t i) = 0;
    virtual uint32_t select1(uint32_t i) = 0;
    virtual uint64_t size_in_bits() = 0;
};

class NaiveBV : public AbstractBV
{
private:
    vector<uint32_t> select0s;
    vector<uint32_t> select1s;

public:
    NaiveBV(vector<bool> &bv)
    {
        for (uint32_t i = 0; i < bv.size(); i++)
        {
            if (bv[i] == 0)
                select0s.push_back(i);
            else
                select1s.push_back(i);
        }
    }

    uint32_t select0(uint32_t i) { return select0s[i - 1]; }
    uint32_t select1(uint32_t i) { return select1s[i - 1]; }

    uint64_t size_in_bits()
    {
        return (select0s.size() + select1s.size()) * 32;
    }
};

class AbstractPredecessor
{
public:
    virtual uint64_t pred(uint64_t x) = 0;
    virtual uint64_t size_in_bits() = 0;
};

class EliasFano : public AbstractPredecessor
{
private:
    AbstractBV *upper_half_bv;
    vector<bool> l_bv;

    uint32_t l_bits;

    uint64_t max_elem;
    uint64_t min_elem;

    uint64_t ith_elem(uint32_t i)
    {
        uint64_t upper_half = upper_half_bv->select1(i + 1) - i;
        uint64_t lower_half = 0;

        for (uint32_t j = (l_bits * i); j < (l_bits * (i + 1)); j++)
        {
            lower_half *= 2;
            lower_half += l_bv[j];
        }

        return (upper_half << l_bits) + lower_half;
    }

public:
    EliasFano(vector<uint64_t> &v)
    {
        min_elem = v.front();
        max_elem = v.back();

        uint32_t u_bits = ceil(log2(v.size()));
        l_bits = ceil(log2(v.back() + 1) - log2(v.size()));

        vector<bool> u_bv = vector<bool>(pow(2, u_bits) - 1 + v.size()); // TODO might need 2n+1, not sure
        l_bv = vector<bool>(v.size() * l_bits);
        for (uint32_t i = 0; i < v.size(); i++)
        {
            uint64_t upper = v[i] >> l_bits;
            u_bv[upper + i] = 1;

            for (uint32_t j = 0; j < l_bits; j++)
                l_bv[(i * l_bits) + j] = (v[i] >> (l_bits - j - 1)) & 1;
        }
        upper_half_bv = new NaiveBV(u_bv);
    }

    uint64_t pred(uint64_t x)
    {
        if (x >= max_elem)
            return max_elem;
        uint64_t x_upper_half = x >> l_bits;

        int32_t p = x_upper_half ? upper_half_bv->select0(x_upper_half) : -1;
        int32_t next_p = upper_half_bv->select0(x_upper_half + 1);

        int32_t left = p - x_upper_half + 1;
        int32_t right = next_p - x_upper_half;

        uint64_t result = left ? ith_elem(left - 1) : UINT64_MAX;
        while (left <= right)
        {
            int32_t mid = left + (right - left) / 2;
            auto elem = ith_elem(mid);
            if (ith_elem(mid) <= x)
            {
                result = elem;
                left = mid + 1;
            }
            else
            {
                right = mid - 1;
            }
        }

        return result;
    }

    uint64_t size_in_bits()
    {
        uint64_t size_in_bits = upper_half_bv->size_in_bits();
        size_in_bits += l_bv.size();
        size_in_bits += 32 + 64 + 64;
        return size_in_bits;
    }
};

enum class RMQ_Algorithm
{
    NAIVE,
    LOGLINEAR,
    LINEAR
};

void print_result(std::string algo, milliseconds time_in_ms, uint64_t space_in_bits)
{
    cout << "RESULT algo=" << algo << " name=atalay_donat time=" << time_in_ms.count() << " space=" << space_in_bits << endl;
}

std::tuple<milliseconds, uint64_t> run_rmq(vector<uint64_t> &v, ifstream &input_file, ofstream &output_file, RMQ_Algorithm rmq_algo = RMQ_Algorithm::LINEAR)
{
    AbstractRMQ *rmq_ds;
    std::string algo_name;

    auto t1 = high_resolution_clock::now();
    switch (rmq_algo)
    {
    case RMQ_Algorithm::NAIVE:
        rmq_ds = new NaiveRMQ(v);
        algo_name = "naivermq";
        break;

    case RMQ_Algorithm::LOGLINEAR:
        rmq_ds = new LogLinearRMQ(v);
        algo_name = "loglinearrmq";
        break;

    case RMQ_Algorithm::LINEAR:
        rmq_ds = new LinearRMQ(v);
        algo_name = "rmq";
    }

    auto t2 = high_resolution_clock::now();
    milliseconds time_in_ms = duration_cast<milliseconds>(t2 - t1);

    uint32_t s, e;
    string line;

    while (input_file >> line)
    {
        size_t comma_idx = line.find(',');
        uint32_t s = stoi(line.substr(0, comma_idx));
        uint32_t e = stoi(line.substr(comma_idx + 1));

        auto t1 = high_resolution_clock::now();
        uint64_t res = rmq_ds->rmq(s, e);
        auto t2 = high_resolution_clock::now();
        time_in_ms += duration_cast<milliseconds>(t2 - t1);
        output_file << res << endl;
    }

    return make_tuple(time_in_ms, rmq_ds->size_in_bits());
}

std::tuple<milliseconds, uint64_t> run_predecessor(vector<uint64_t> &v, ifstream &input_file, ofstream &output_file)
{
    auto t1 = high_resolution_clock::now();
    AbstractPredecessor *pd_ds = new EliasFano(v);
    auto t2 = high_resolution_clock::now();
    milliseconds time_in_ms = duration_cast<milliseconds>(t2 - t1);

    uint64_t x;
    while (input_file >> x)
    {
        auto t1 = high_resolution_clock::now();
        uint64_t res = pd_ds->pred(x);
        auto t2 = high_resolution_clock::now();
        time_in_ms += duration_cast<milliseconds>(t2 - t1);
        output_file << res << endl;
    }

    return make_tuple(time_in_ms, pd_ds->size_in_bits());
}

int main(int argc, char *argv[])
{
    ifstream input_file(argv[2]);
    if (!input_file.is_open())
        throw runtime_error("Unable to open input file "s + argv[2]);
    ofstream output_file(argv[3]);
    if (!output_file.is_open())
        throw runtime_error("Unable to open output file "s + argv[3]);

    uint32_t n;
    input_file >> n;

    vector<uint64_t> v(n);
    for (uint32_t i = 0; i < n; i++)
        input_file >> v[i];

    milliseconds time_in_ms;
    uint64_t space_in_bits;

    if (argv[1] == "pd"s)
        tie(time_in_ms, space_in_bits) = run_predecessor(v, input_file, output_file);
    else if (argv[1] == "rmq"s)
        tie(time_in_ms, space_in_bits) = run_rmq(v, input_file, output_file);
    else if (argv[1] == "loglinearrmq"s)
        tie(time_in_ms, space_in_bits) = run_rmq(v, input_file, output_file, RMQ_Algorithm::LOGLINEAR);
    else if (argv[1] == "naivermq"s)
        tie(time_in_ms, space_in_bits) = run_rmq(v, input_file, output_file, RMQ_Algorithm::NAIVE);
    else
        throw invalid_argument("The first parameter must either be 'pd', 'rmq', 'loglinearrmq' or 'naivermq'.");

    print_result(argv[1], time_in_ms, space_in_bits);
}