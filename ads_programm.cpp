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
 * @brief Abtract class for data structures which support Range Minimum Queries.
 */
class AbstractRMQ
{
public:
    /**
     * @brief Returns the RMQ in the given range.
     *
     * @param s The start of the range
     * @param e The end of the range (inclusive)
     */
    virtual uint32_t rmq(uint32_t s, uint32_t e) = 0;
    /**
     * @brief Returns the cumulative size of all the instance variables in bits.
     */
    virtual uint64_t size_in_bits() = 0;
};

/**
 * @brief Range Minimum Query data structure which stores the result of all the possible ranges, which results in O(n^2) space.
 */
class NaiveRMQ : public AbstractRMQ
{
private:
    /**
     * @brief RMQ solutions for all possible ranges, where RMQ(s, e) = rmq_solutions[s][e - s - 1].
     */
    vector<vector<uint32_t>> rmq_solutions;

public:
    /**
     * @brief Construct a new NaiveRMQ instance with the given vector.
     *
     * @param v Vector on which RMQ values will be queried.
     */
    NaiveRMQ(vector<uint64_t> &v)
    {
        uint32_t n = v.size();

        /*
        For a given query RMQ(s, e), s must be <= e, so space is not allocated for s < e in rmq_soluions.
        Furthermore, RMQ(s, s) = s, so this is also not needed to be stored. Overall, each s needs to store
        n-s-1 solutions. Therefore, RMQ(s, e) = rmq_solutions[s][e - s - 1]
        */

        // the last element in the vector does not store any solution, so the size is n-1.
        rmq_solutions = vector<vector<uint32_t>>(n - 1);

        // iteratively compute the RMQ solutions and store them
        for (uint32_t i = 0; i < n - 1; i++)
        {
            rmq_solutions[i] = vector<uint32_t>(n - i - 1);

            uint32_t cur_min_idx = i;
            for (uint32_t j = i + 1; j < n; j++)
            {
                // invariant: cur_min_idx = RMQ(i, j - 1)
                // RMQ(i, j) = j if v[j] < v[RMQ(i, j - 1)], else RMQ(i, j - 1)
                if (v[j] < v[cur_min_idx])
                    cur_min_idx = j;
                rmq_solutions[i][j - i - 1] = cur_min_idx;
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
        // only rmq_solutions is stored
        uint64_t size_in_bits = 0;
        for (auto &rmq_sol : rmq_solutions)
            size_in_bits += rmq_sol.size() * 32;

        return size_in_bits;
    }
};

/**
 * @brief Range Minimum Query data structure which stores the result of ranges whose length is a power of 2, which results in O(nlogn) space.
 */
class LogLinearRMQ : public AbstractRMQ
{
private:
    /**
     * @brief RMQ solutions for all possible ranges whose length is a power of 2 and greater than 1, where RMQ(s, s + 2^l - 1) = rmq_solutions[s][l - 1].
     */
    vector<vector<uint32_t>> rmq_solutions;
    /**
     * @brief Pointer to the underlying vector
     */
    vector<uint64_t> *v;

public:
    LogLinearRMQ(vector<uint64_t> &v)
    {
        uint32_t n = v.size();

        // RMQ(s, s) = RMQ(s, s + 2^0 - 1) = s, so solutions for l=0 will not be stored. Therefore, RMQ(s, s + 2^l - 1) = rmq_solutions[s][l - 1]

        // the last element in the vector does not store any solution, so the size is n-1.
        rmq_solutions = vector<vector<uint32_t>>(n - 1);
        this->v = &v;

        // dynamic programming is used to construct rmq_solutions as described in the lecture in O(nlogn) time
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

        // Answers RMQ using two subqueries as described in the lecture which can be extracted from rmq_solutions.

        /*
        log2(e-s-1) rounded down as described in the lecture does not work because it might be the case that the ranges do not overlap.
        Example: s=0, e=2: log(e-s-1)=log2(1)=0, so m1 = RMQ(s, s+2^0-1) = RMQ(s, s) = RMQ(0, 0) and m2 = RMQ(e-2^0+1, e) = RMQ(e, e) = RMQ(2, 2).
        The ranges do not overlap.
        */
        uint32_t l = log2(e - s);

        // m1 = RMQ(s, s + 2^l - 1)
        uint32_t m1 = (l > 0) ? rmq_solutions[s][l - 1] : s;
        // m2 = RMQ(e - 2^l + 1, e)
        uint32_t m2 = (l > 0) ? rmq_solutions[e - pow(2, l) + 1][l - 1] : e - pow(2, l) + 1;

        return ((*v)[m1] < (*v)[m2]) ? m1 : m2;
    }

    uint64_t size_in_bits()
    {
        // the vector needs to be stored for comparing the answers of the two subqueries
        uint64_t size_in_bits = (*v).size() * 64;

        // size of rmq_solutions
        for (auto &rmq_sol : rmq_solutions)
            size_in_bits += rmq_sol.size() * 32;

        return size_in_bits;
    }
};

/**
 * @brief Range Minimum Query data structure which stores RMQ of possible Cartesian trees, which results in O(n) space.
 */
class LinearRMQ : public AbstractRMQ
{
private:
    uint32_t block_size;
    vector<uint64_t> *v;

    /**
     * @brief The minimum value within each block
     */
    vector<uint64_t> min_within_block;

    /**
     * @brief The index of the minimum value within each block
     */
    vector<uint32_t> min_idx_within_block;
    /**
     * @brief RMQ data structure for spanning block queries
     */
    AbstractRMQ *query_spanning_block_rmq_ds;

    /**
     * @brief Cartesian trees of each block in Cartesian tree signature representation
     */
    vector<vector<bool>> c_trees;
    /**
     * @brief RMQ solutions of each Cartesian tree at each possible range. For a given cartesian tree c_tree, RMQ(s, e) = c_tree_start_end_rmqs[c_tree][s][e - 1]
     */
    unordered_map<vector<bool>, vector<vector<uint32_t>>> c_tree_start_end_rmqs;

    /**
     * @brief Construct the Cartesian tree of the block in the vector in the given range in the Cartesian tree signature representation
     * described in the first paragraph of Sec. 3.7 in https://drops.dagstuhl.de/opus/volltexte/2019/10487/pdf/LIPIcs-CPM-2019-16.pdf (without the
     * redundant 0 in the first bit of the representation).
     *
     * @param start_idx start index of the range
     * @param end_idx end index of the range (exclusive)
     *
     * @return Cartesian tree signature representation of the Cartesian tree of the block
     */
    vector<bool> construct_c_tree(vector<uint64_t> &vec, uint32_t start_idx, uint32_t end_idx)
    {
        vector<bool> repr;
        vector<uint64_t> stack;

        /*
        The Cartesian tree signature is constructed iteratively as follows:
        The rightmost path in the Cartesian tree (computed until now) is stored in the stack, where the leaf is at the top.
        Pop the stack until the current value is not smaller than the top of the stack. In the representation, append 1 for each time the stack is popped, and 0 at the end.
        Then, push the current value to the stack. The redundant 0 at the start of the representation is removed in this representation.
        */

        stack.push_back(vec[start_idx]);

        for (size_t i = start_idx + 1; i < end_idx; i++)
        {
            while ((!stack.empty()) && (vec[i] < stack.back()))
            {
                repr.push_back(1);
                stack.pop_back();
            }

            repr.push_back(0);
            stack.push_back(vec[i]);
        }

        /*
        When running LinearRMQ, it might be the case that the range size is less than block size (at the last block). This causes problems when indexing with this tree, so
        0s are appended to the representation. This represents the Cartesian tree of a block with the required block size where the RMQ solutions are the same as the current
        block.
        */
        for (uint32_t i = end_idx - start_idx; i < block_size; i++)
            repr.push_back(0);

        return repr;
    }

    /**
     * @brief Store RMQ solutions for each possible Cartesian tree with the block size and for each possible range
     */
    void construct_c_tree_start_end_rmqs()
    {
        vector<uint64_t> vec(block_size);

        // initialize vector with values from 0 to n-1
        for (uint32_t i = 0; i < block_size; i++)
            vec[i] = i;

        /*
        For each permutation of the numbers 0..block size - 1, calculate cartesian tree and RMQ for all possible ranges with NaiveRMQ.
        Block sizes tend to be very small, so it is not worth to use LoglinearRMQ here.
        */
        do
        {
            AbstractRMQ *rmq_ds = new NaiveRMQ(vec);
            vector<bool> c_tree = construct_c_tree(vec, 0, block_size);

            // pass if RMQ values for this cartesian tree are already calculated
            if (c_tree_start_end_rmqs.find(c_tree) != c_tree_start_end_rmqs.end())
                continue;

            /*
            For a given query RMQ(s, e), s must be <= e, so space is not allocated for s < e in rmq_soluions.
            Furthermore, RMQ(s, s) = s, so this is also not needed to be stored. Overall, each s needs to store
            n-s-1 solutions. Therefore, for a given cartesian tree c_tree: RMQ(s, e) = c_tree_start_end_rmqs[c_tree][s][e - s - 1]
            */

            c_tree_start_end_rmqs[c_tree] = vector<vector<uint32_t>>(block_size - 1);

            for (uint32_t i = 0; i < block_size - 1; i++)
            {
                c_tree_start_end_rmqs[c_tree][i] = vector<uint32_t>(block_size - i);
                for (uint32_t j = i + 1; j < block_size; j++)
                    c_tree_start_end_rmqs[c_tree][i][j - i - 1] = rmq_ds->rmq(i, j);
            }
        } while (next_permutation(vec.begin(), vec.end()));
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

        construct_c_tree_start_end_rmqs();

        // for each block, store the smallest value, its index and its cartesian tree
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

            vector<bool> block_c_tree = construct_c_tree(v, start_idx, end_idx);
            c_trees[i] = block_c_tree;
        }

        // recursively create LinearRMQs for the blocks. In the base case, create LogLinearRMQs instead
        // this is drastically superior to using directly the LogLinearRMQ in terms of both time and space
        if (number_of_blocks > 100)
            query_spanning_block_rmq_ds = new LinearRMQ(min_within_block);
        else
            query_spanning_block_rmq_ds = new LogLinearRMQ(min_within_block);
    }

    /**
     * @brief Returns the minimum value between the given blocks
     *
     * @param s_block start block
     * @param e_block end block (inclusive)
     */
    uint32_t spanning_block_rmq(uint32_t s_block, uint32_t e_block)
    {
        uint32_t min_block_idx = query_spanning_block_rmq_ds->rmq(s_block, e_block);
        return min_idx_within_block[min_block_idx];
    }

    /**
     * @brief Returns the minimum value within the given block and the given range
     *
     * @param block_idx index of the block
     * @param s start of the range
     * @param e end of the range (inclusive)
     */
    uint32_t within_block_rmq(uint32_t block_idx, uint32_t s, uint32_t e)
    {
        if (s == e)
            return (block_idx * block_size) + s;
        return (block_idx * block_size) + c_tree_start_end_rmqs[c_trees[block_idx]][s][e - s - 1];
    }

    uint32_t rmq(uint32_t s, uint32_t e)
    {
        // the blocks of the start and end index
        uint32_t s_block = (double)s / block_size;
        uint32_t e_block = (double)e / block_size;

        // in each case, return the index which corresponds to the minimum value out of other indices from the other subqueries
        if (s_block == e_block)
        {
            // the range is in one block, just use RMQ within block
            return within_block_rmq(s_block, s % block_size, e % block_size);
        }
        else if ((s % block_size == 0) && ((e % block_size == block_size - 1) || (e == v->size() - 1)))
        {
            // s is the first index in s_block and e is the last index in e_block, just use spanning block RMQ
            return spanning_block_rmq(s_block, e_block);
        }
        else if (s % block_size == 0)
        {
            // s is first index in s_block but e is not the last index in e_block, so spanning block RMQ from s_block to e_block -1 + within block RMQ in e_block
            uint32_t within_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block, e_block - 1);
            return ((*v)[within_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_block_min_idx : spanning_block_min_idx;
        }
        else if ((e % block_size == block_size - 1) || (e == v->size() - 1))
        {
            // s is not first index in s_block but e is the last index in e_block, so within block RMQ in s_block + spanning block RMQ from s_block + 1 to e_block
            uint32_t within_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t spanning_block_min_idx = spanning_block_rmq(s_block + 1, e_block);
            return ((*v)[within_block_min_idx] < (*v)[spanning_block_min_idx]) ? within_block_min_idx : spanning_block_min_idx;
        }
        else if (s_block + 1 == e_block)
        {
            // s is not first index in s_block, e is not last index in e_block and s_block is just before e_block, so just use within block RMQ in both blocks
            uint32_t within_s_block_min_idx = within_block_rmq(s_block, s % block_size, block_size - 1);
            uint32_t within_e_block_min_idx = within_block_rmq(e_block, 0, e % block_size);
            return ((*v)[within_s_block_min_idx] < (*v)[within_e_block_min_idx]) ? within_s_block_min_idx : within_e_block_min_idx;
        }
        else
        {
            // within block RMQ in s_block + spanning block RMQ from s_block + 1 to e_block - 1 + within block RMQ in e_block
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

/**
 * @brief Abstract bitvector data structure supporting select queries
 */
class AbstractBV
{
public:
    /**
     * @brief Returns the index of the i-th 0 in the bitvector.
     */
    virtual uint32_t select0(uint32_t i) = 0;
    /**
     * @brief Returns the index of the i-th 1 in the bitvector.
     */
    virtual uint32_t select1(uint32_t i) = 0;
    /**
     * @brief Returns the size of the data structure instance in bits.
     */
    virtual uint64_t size_in_bits() = 0;
};

/**
 * @brief Concrete implementation of a bitvector data structure supporting select queries
 */
class ConcreteBV : public AbstractBV
{
private:
    vector<uint32_t> select0s;
    vector<uint32_t> select1s;

public:
    ConcreteBV(vector<bool> &bv)
    {
        for (uint32_t i = 0; i < bv.size(); i++)
        {
            if (bv[i] == 0)
                select0s.push_back(i);
            else
                select1s.push_back(i);
        }
    }

    // subtract by 1 because of 0-indexed vectors
    uint32_t select0(uint32_t i) { return select0s[i - 1]; }
    uint32_t select1(uint32_t i) { return select1s[i - 1]; }

    uint64_t size_in_bits()
    {
        return (select0s.size() + select1s.size()) * 32;
    }
};

/**
 * @brief Abstract data structure supporting predecessor queries in a vector
 */
class AbstractPredecessor
{
public:
    /**
     * @brief Return the predecessor of the given value
     */
    virtual uint64_t pred(uint64_t x) = 0;
    /**
     * @brief Returns the size of the instance in bits
     */
    virtual uint64_t size_in_bits() = 0;
};

/**
 * @brief Elias Fano algorithm supporting predecessor queries in a vector
 */
class EliasFano : public AbstractPredecessor
{
private:
    /**
     * @brief Bitvector data structure of the upper half supporting select queries
     */
    AbstractBV *upper_half_bv;

    uint32_t u_bv_size;

    /**
     * @brief Bitvector of the lower half
     */
    vector<bool> l_bv;

    /**
     * @brief Number of lower bits of each element 
     */
    uint32_t l_bits;

    /**
     * @brief Number of upper bits of each element 
     */
    uint32_t u_bits;

    /**
     * @brief Returns the i-th element in the vector
     */
    uint64_t ith_elem(uint32_t i)
    {
        // this holds because of how we insert to the upper half bitvector, i+1 because i is 0-indexed
        uint64_t upper_half = upper_half_bv->select1(i + 1) - i;

        // get lower half from the bit representation in l_bv
        uint64_t lower_half = 0;
        for (uint32_t j = (l_bits * i); j < (l_bits * (i + 1)); j++)
        {
            lower_half <<= 1;
            lower_half += l_bv[j];
        }

        return (upper_half << l_bits) + lower_half;
    }

    /**
     * @brief Get the upper half of the given value
     */
    uint64_t get_upper_half(uint64_t x) {
        return (x >> l_bits) & ((((uint64_t) 1) << u_bits) - 1);
    }

public:
    EliasFano(vector<uint64_t> &v)
    {
        // number of bits in upper half and lower half
        u_bits = ceil(log2(v.size()));
        l_bits = ceil(log2(v.back() + 1) - log2(v.size()));

        /*
        The length of upper half bitvector is 2n + 1 in the lecture, but also in the example in the lecture (lecture 4 slide 6), but suppose
        the upper half is 1111 and the position is 9. In the bitvector, this corresponds to a 1 in position 15 + 9 = 24, which does not fit into
        a bitvector of size 2 * 10 + 1 = 21. The problem is that the number of bits upper half uses is rounded up, so the maximum upper half can 
        be 2^(logn rounded up) - 1 = 2^u_bits - 1. The maximum position in the vector can be n - 1, so the maximum position in the bitvector can
        be 2^u_bits - 1 + n - 1. Therefore the size of the upper half bitvector is 2^u_bits - 1 + n.
        */ 
        vector<bool> u_bv = vector<bool>(pow(2, u_bits) - 1 + v.size());
        u_bv_size = u_bv.size();
        // lower bit vector stores l_bits per element in the vector
        l_bv = vector<bool>(v.size() * l_bits);
        for (uint32_t i = 0; i < v.size(); i++)
        {
            uint64_t upper = get_upper_half(v[i]);
            u_bv[upper + i] = 1;

            // store the lower bits one by one in the bitvector l_bv
            for (uint32_t j = 0; j < l_bits; j++)
                l_bv[(i * l_bits) + j] = (v[i] >> (l_bits - j - 1)) & 1;
        }

        upper_half_bv = new ConcreteBV(u_bv);
    }

    uint64_t pred(uint64_t x)
    {
        uint64_t x_upper_half = get_upper_half(x);

        // index of the zeroes in the upper half bitvector between which the predecessor might be in (p is -1 if the upper half is 0, it works with the variables left and right)
        int32_t p = x_upper_half ? upper_half_bv->select0(x_upper_half) : -1;
        int32_t next_p = (x_upper_half != (pow(2, u_bits))) ? upper_half_bv->select0(x_upper_half + 1) : u_bv_size;

        // indices of the search range in the vector
        int32_t left = p - x_upper_half + 1;
        int32_t right = next_p - x_upper_half;

        // if the predecessor is not in the search range, then the predecessor is the element just before the range, or maximum 64-bit integer if there is no predecessor
        uint64_t result = left ? ith_elem(left - 1) : UINT64_MAX;

        // binary search to find the predecessor
        while (left <= right)
        {
            int32_t mid = left + (right - left) / 2;
            uint64_t elem = ith_elem(mid);
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
        size_in_bits += 32 + 32;
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