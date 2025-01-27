#include <iostream>
#include <cstdint>
#include <random>
#include <algorithm>
#include <set>
#include <chrono>
#include <vector>
#include <utility>

#include "parlay/random.h"
#include "parlay/sequence.h"
#include "knuth_shuffle.h"
#include "parlay/primitives.h"


void Gen1 (uint64_t n) {
    parlay::random_generator gen;
    std::uniform_int_distribution<uint64_t> dis(0,1);
    auto seq = parlay::tabulate(n, [&] (uint64_t i) {
        return i;
    });

    std::random_device rd;
    uint64_t salt = rd();

    auto random_bits = parlay::tabulate(n, [&] (auto i) {
        auto r = gen[i ^ salt];
        return dis(r);
    });

    auto S = parlay::scan(random_bits);
    auto t = S.second;
    auto C = parlay::sequence<uint64_t>(n);

    parlay::parallel_for(0, n, [&] (auto i) {
        if (random_bits[i] == 1) {
            C[S.first[i]] = seq[i];
        } else {
            C[t + i - S.first[i]] = seq[i];
        }
    });
}

void Gen2 (uint64_t n) {
    std::random_device rd;
    uint64_t salt = rd();
    size_t mid = n / 2;
    auto seq = parlay::random_permutation(n, salt);
    auto A = parlay::make_slice(seq.begin(), seq.begin() + mid);
    auto B = parlay::make_slice(seq.begin() + mid, seq.end());

    auto C = parlay::integer_sort(A);
    auto D = parlay::integer_sort(B);


    // size_t mid = n / 2;

    // parlay::parallel_do(
    //     [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin(), seq.begin() + mid)); },
    //     [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin() + mid, seq.end())); }
    // );

    // return std::move(seq);
}

// parlay::sequence<uint64_t> Gen3 (uint64_t n) {
//     auto seq = parlay::tabulate(n, [&] (auto i) {
//         return i;
//     });

//     random_shuffle(seq);
//     size_t mid = n / 2;

//     parlay::parallel_do(
//         [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin(), seq.begin() + mid)); },
//         [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin() + mid, seq.end())); }
//     );

//     return std::move(seq);
// }


int main(int argc, char *argv[]) {

   auto n = static_cast<uint64_t>(atoi(argv[1]));
   auto start_1 = std::chrono::high_resolution_clock::now();
   Gen1(n);
   auto end_1 = std::chrono::high_resolution_clock::now();
   auto time_1 = std::chrono::duration_cast<std::chrono::microseconds>(end_1-start_1);

   auto start_2 = std::chrono::high_resolution_clock::now();
   Gen2(n);
   auto end_2 = std::chrono::high_resolution_clock::now();
   auto time_2 = std::chrono::duration_cast<std::chrono::microseconds>(end_2-start_2);

   
//    auto start_3 = std::chrono::high_resolution_clock::now();
//    auto R_3 = Gen3(n);
//    auto end_3 = std::chrono::high_resolution_clock::now();
//    auto time_3 = std::chrono::duration_cast<std::chrono::microseconds>(end_3-start_3);


   std::cout << "\n\nGen1 time in microseconds: " << time_1.count() << "\nGen2 time in microseconds: " << time_2.count() <<  "\n\n";
   return 0; 
}