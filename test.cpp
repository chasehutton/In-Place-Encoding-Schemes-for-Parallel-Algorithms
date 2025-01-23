#include <iostream>
#include <cstdint>
#include <random>
#include <algorithm>
#include <set>
#include <chrono>

#include "utils.h"
#include "merge.h"



#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"

// This file contains GPT generated test code

parlay::sequence<uint32_t> generateUniqueTwoSortedHalves(uint32_t size) {
    if (size % 2 != 0) {
        throw std::invalid_argument("Size must be even.");
    }
    size_t half = size / 2;

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dis(1, 3000000000);

    std::set<uint32_t> globalSet;
    while (globalSet.size() < size) {
        globalSet.insert(dis(gen));
    }

    parlay::sequence<uint32_t> allVals(globalSet.begin(), globalSet.end());
    std::shuffle(allVals.begin(), allVals.end(), gen);

    parlay::sequence<uint32_t> result(size);
    std::copy(allVals.begin(), allVals.begin() + half, result.begin());
    std::copy(allVals.begin() + half, allVals.end(),result.begin() + half);

    std::sort(result.begin(),  result.begin() + half);
    std::sort(result.begin() + half,  result.end());

    return result;
}

int main() {
    size_t size = 20480000;

    parlay::sequence<uint32_t> testSequence = generateUniqueTwoSortedHalves(size);
    std::cout << "Testing...\n\n";
     // measure memory before

    auto start = std::chrono::high_resolution_clock::now();

    Merge(testSequence, 2048);

    auto end = std::chrono::high_resolution_clock::now();

    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);

    std::cout << "Time taken by IN-PLACE Merge: "
         << time.count() << " microseconds" << std::endl;


    auto A = testSequence.subseq(0, testSequence.size()/2 - 1);
    auto B = testSequence.subseq(testSequence.size()/2, testSequence.size() - 1);

    start = std::chrono::high_resolution_clock::now();

    auto r = parlay::merge(A, B);

    end = std::chrono::high_resolution_clock::now();

    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);

    std::cout << "Time taken by Merge: "
         << time.count() << " microseconds" << std::endl;


    return 0;
}