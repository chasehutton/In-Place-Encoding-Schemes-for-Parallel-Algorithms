#pragma once

#include <cstdint>


#include "parlay/sequence.h"
#include "parlay/parallel.h"

/////////////////////////////////////////////////////////////////////
// Basic block utility functions each overloaded to work with      //
// Parlay Slices.                                                  //
/////////////////////////////////////////////////////////////////////


inline void SwapBlock(parlay::sequence<uint32_t>& seq,
                uint32_t block1_start, uint32_t block1_end, uint32_t block2_start, uint32_t block2_end) {
    uint32_t size = block1_end - block1_start + 1;
    // assert(size == block2_end - block2_start + 1);
    // assert(size % 2 == 0);

    uint32_t t = 0;
    for (int i = 0; i < size; i++) {
        t = seq[block1_start + i];
        seq[block1_start + i] = seq[block2_start + i];
        seq[block2_start + i] = t;
    }
}

inline void SwapBlock(parlay::slice<uint32_t*, uint32_t*> block1, 
                parlay::slice<uint32_t*, uint32_t*> block2) {
    uint32_t size = block1.size();
    // assert(size == block2.size());

    parlay::parallel_for(0, size, [&] (auto i) {
        std::swap(block1[i], block2[i]);
    });
}

inline uint32_t ReadBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end) {
    uint32_t size = end - start + 1;
    // assert(size % 2 == 0);
    // assert(size/2 <= 32);

    uint32_t result = 0;
    uint8_t bit_pos = 0;
    for(int i = 0; i < (size/2); i++) {
        if (block[start + 2*i] > block[start + 2*i + 1]) result |= (1U << bit_pos);
        bit_pos++;
    }

    return result;
}


inline void WriteBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end, uint32_t value) {
    uint32_t size = end - start + 1;
    // assert(size % 2 == 0);
    // assert(size/2 <= 32);
    for (int i = 0; i < (size/2); i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = block[start + 2*i];
        uint32_t& second = block[start + 2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}

inline void WriteBlock(parlay::slice<uint32_t*, uint32_t*> block, uint32_t value) {
    uint32_t size = block.size();
    // assert(size % 2 == 0);
    // assert(size/2 <= 32);
    for (int i = 0; i < (size/2); i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = block[2*i];
        uint32_t& second = block[2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}


