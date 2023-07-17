# Advanced Data Structure Project

This project contains succinct implementations for Range Minimum Query and Predecessor queries on a vector containing 64-bit elements (and the length of the vector should be 32-bits) as follows:

- Elias Fano Predecessor
- Naive RMQ (O(n^2) space)
- Loglinear space RMQ
- Linear space RMQ with Cartesian trees

## Compiling 

```
g++ -std=c++20 -O3 ads_programm.cpp -o ads_programm
```

## Running

The input file of each command should correspond to the style described in the project description.

- Predecessor
```
./ads_programm pd <input_file> <output_file>
```

- Linear RMQ
```
./ads_programm rmq <input_file> <output_file>
```

- Loglinear RMQ
```
./ads_programm loglinearrmq <input_file> <output_file>
```

- Naive RMQ
```
./ads_programm naivermq <input_file> <output_file>
```