#pragma once

#include <cstdint>


#include "parlay/sequence.h"
#include "parlay/parallel.h"

void SwapBlockN(parlay::sequence<uint64_t>& block1, parlay::sequence<uint64_t>& block2, uint64_t granularity) {
    uint64_t size = block1.size();
    assert(size == block2.size());
    assert(size % 2 == 0);

    parlay::parallel_for(0, size, [&] (auto i) {
        std::swap(block1[i], block2[i]);
    }, granularity);
}

uint32_t ReadBlock64(const parlay::sequence<uint64_t>& block, uint64_t start, uint64_t end) {
    assert((end - start) == 64);
    uint32_t result = 0;

    for(int i = 0; i < 32; i++) {
        result <<= 1;
        result |= (block[start + 2*i] < block[start + 2*i + 1]) ? 0 : 1;
    }

    return result;
}

void WriteBlock64(parlay::sequence<uint64_t>& block, uint64_t start, uint64_t end, const uint32_t value) {
    assert((end - start) == 64);
    for (int i = 0; i < 32; i++) {
        bool b = (value << i) & 1U;
        uint64_t& first = block[start + 2*i];
        uint64_t& second = block[start + 2*i + 1];
        if ((b == 0 && first > second) || (b == 1 && first < second)) {
            std::swap(first, second);
        } 
    }
}


