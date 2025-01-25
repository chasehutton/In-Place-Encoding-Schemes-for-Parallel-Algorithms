#pragma once

#include <cstdint>
#include <random>
#include <iostream>
#include <chrono>


#include "utils.h"
#include "buffer.h"

#include "parlay/sequence.h"
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"

void set_up(parlay::sequence<uint32_t>& seq, uint32_t c) {
	buffer 
}


void Merge(parlay::sequence<uint32_t>& seq, uint32_t b) {
    assert(seq.size()/2 % b == 0);
    assert(seq.size() > b);
    assert(b % 2 == 0);
    assert(b >= 5*SEGMENT_SIZE);

    bool* flag = (bool*) std::malloc(sizeof(bool));
    *flag = false;
    std::cout << "Setting Up...\n\n";
    auto start = std::chrono::high_resolution_clock().now();
    SetUp(seq, b);
    auto end = std::chrono::high_resolution_clock().now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "Time for SetUp: " << time.count() << " \n";
    std::cout << "End Sorting...\n\n";
    start = std::chrono::high_resolution_clock().now();
    EndSort(seq, b, flag);
    end = std::chrono::high_resolution_clock().now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "Time for EndSort: " << time.count() << " \n";
    std::cout << "Seq Sorting...\n\n";
    start = std::chrono::high_resolution_clock().now();
    SeqSort(seq, 0, (uint32_t)seq.size(), b);
    end = std::chrono::high_resolution_clock().now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "Time for SeqSort: " << time.count() << " \n";

    std::free(flag);
}

