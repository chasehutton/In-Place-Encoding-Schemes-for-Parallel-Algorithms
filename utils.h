#pragma once

#include <cstdint>


#include "parlay/sequence.h"
#include "parlay/parallel.h"

void SwapBlockN(parlay::sequence<uint32_t>& block1, parlay::sequence<uint32_t>& block2, uint32_t granularity) {
    uint32_t size = block1.size();
    assert(size == block2.size());
    assert(size % 2 == 0);

    parlay::parallel_for(0, size, [&] (auto i) {
        std::swap(block1[i], block2[i]);
    }, granularity);
}

uint32_t ReadBlock32(const parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end) {
    assert((end - start) == 64);
    uint32_t result = 0;

    for(int i = 0; i < 32; i++) {
        result <<= 1;
        result |= (block[start + 2*i] < block[start + 2*i + 1]) ? 0 : 1;
    }

    return result;
}

void WriteBlock32(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end, const uint32_t value) {
    assert((end - start) == 64);
    for (int i = 0; i < 32; i++) {
        uint32_t b = value & (1U << i);
        uint32_t& first = block[start + 2*i];
        uint32_t& second = block[start + 2*i + 1];
        if ((b == 0 && first > second) || (b == 1 && first < second)) {
            std::swap(first, second);
        } 
    }
}


