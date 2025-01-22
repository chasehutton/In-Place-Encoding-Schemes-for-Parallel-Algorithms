#pragma once

#include <cstdint>


#include "parlay/sequence.h"
#include "parlay/parallel.h"

/////////////////////////////////////////////////////////////////////
// Basic block utility functions each overloaded to work with      //
// Parlay Slices.                                                  //
/////////////////////////////////////////////////////////////////////


inline void SwapBlock(parlay::sequence<uint32_t>& block1, parlay::sequence<uint32_t>& block2,
                uint32_t block1_start, uint32_t block1_end, uint32_t block2_start, uint32_t block2_end, uint32_t granularity) {
    uint32_t size = block1_end - block1_start;
    assert(size == block2_end - block2_start);

    parlay::parallel_for(0, size, [&] (auto i) {
        std::swap(block1[block1_start + i], block2[block2_start + i]);
    }, granularity);
}

inline void SwapBlock(parlay::slice<uint32_t*, uint32_t*> block1, 
                parlay::slice<uint32_t*, uint32_t*> block2,
                uint32_t granularity) {
    uint32_t size = block1.size();
    assert(size == block2.size());

    parlay::parallel_for(0, size, [&] (auto i) {
        std::swap(block1[i], block2[i]);
    }, granularity);
}

inline uint32_t ReadBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end) {
    uint32_t size = end - start;
    assert(size % 2 == 0);
    uint32_t result = 0;
    for(int i = 0; i < (size/2); i++) {
        result <<= 1;
        result |= (block[start + 2*i] < block[start + 2*i + 1]) ? 0 : 1;
    }

    return result;
}

inline uint32_t ReadBlock(parlay::slice<uint32_t*, uint32_t*> block) {
    uint32_t size = block.size();
    assert(size % 2 == 0);
    uint32_t result = 0;
    for(int i = 0; i < (size/2); i++) {
        result <<= 1;
        result |= (block[2*i] < block[2*i + 1]) ? 0 : 1;
    }

    return result;
}

inline void WriteBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end, uint32_t value) {
    uint32_t size = end - start;
    assert(size % 2 == 0);
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
    assert(size % 2 == 0);
    for (int i = 0; i < (size/2); i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = block[2*i];
        uint32_t& second = block[2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}


