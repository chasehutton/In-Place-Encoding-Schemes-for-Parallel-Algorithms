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

#define SEGMENT_SIZE 64

inline uint32_t log2(uint32_t x) {
	uint32_t e = -1;
 	for(int i = 0; i < 32; i++) {
 		if (x & 1) e = i;
 		x >> 1;
 	}
 	return e;
}

void set_up(parlay::sequence<uint32_t>& seq, buffer& Buffer, uint32_t c) {
	
}


void Merge(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B,  uint32_t b) {
    assert(seq.size()/2 % b == 0);
    assert(seq.size() > b);
    assert(b % 2 == 0);
    assert(b >= 5*SEGMENT_SIZE);

	uint32_t n = seq.size();
    uint32_t r = n/(segment*log2(n));

    r + r log2(n)

    buffer Buffer = new buffer(parlay::make_slice());

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

