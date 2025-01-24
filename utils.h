#pragma once

#include <cstdint>


#include "parlay/sequence.h"
#include "parlay/parallel.h"

///////////////////////////////////////////////////////////
// Blocks are given by seq[block_start, block_end).      //
///////////////////////////////////////////////////////////

// Assumes block1_end - block1_start = block2_end - block2_start
// Assumes size is divisble by 2
inline void SwapBlock(parlay::sequence<uint32_t>& seq, uint32_t block1_start, uint32_t block1_end,
                      uint32_t block2_start, uint32_t block2_end) {

    uint32_t size = block1_end - block1_start;
    for (int i = 0; i < size; i++) {
        std::swap(seq[block1_start + i], seq[block2_start + i]);
    }
}

inline void SwapBlock(parlay::slice<uint32_t*, uint32_t*> block1, 
                parlay::slice<uint32_t*, uint32_t*> block2) {
    
    uint32_t size = block1.size();
    for (int i = 0; i < size; i++) {
        std::swap(block1[i], block2[i]);
    }
}

// Assumes size is divisible by 2 and less than or equal to 64
inline uint32_t ReadBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end) {
    uint32_t size = end - start;
    uint32_t result = 0;
    uint8_t bit_pos = 0;
    for(int i = 0; i < size/2; i++) {
        if (block[start + 2*i] > block[start + 2*i + 1]) result |= (1U << bit_pos);
        bit_pos++;
    }

    return result;
}

inline uint32_t ReadBlock(parlay::slice<uint32_t*, uint32_t*> block) {
    uint32_t size = block.size();
    uint32_t result = 0;
    uint32_t bit_pos = 0;
    for (int i = 0; i < size/2; i++) {
        if (block[2*i] > block[2*i + 1]) result |= (1U << bit_pos);
        bit_pos++;
    }

    return result;
}

// Assumes size is divisible by 2 and less than or equal to 64
inline void WriteBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end, uint32_t value) {
    uint32_t size = end - start;
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
    for (int i = 0; i < (size/2); i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = block[2*i];
        uint32_t& second = block[2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}

// Simple BubbleSort
inline void BubbleSort(parlay::slice<uint32_t*, uint32_t*> A, parlay::slice<uint32_t*, uint32_t*> B) {
    auto nA = A.size();
    auto nB = B.size();
    auto n  = nA + nB;
    bool swapped;
    for (int i = 0; i < n; i++) {
        swapped = false;
        for (int j = 0; j < n-i-1; j++) {
            if (j+1 < nA) {
                if (A[j] > A[j+1]) {
                    std::swap(A[j], A[j+1]);
                    swapped = true;
                }
            } else if (j >= nA) {
                if (B[j - nA] > B[j+1 - nA]) {
                    std::swap(B[j - nA], B[j+1 - nA]);
                    swapped = true;
                }
            } else {
                if (A[nA - 1] > B[0]) {
                    std::swap(A[nA - 1], B[0]);
                    swapped = true;
                }
            }
        }

        if (!swapped) return;
    }
}

// Sorting seq[start, end)
inline void BubbleSort(parlay::sequence<uint32_t>& seq, uint32_t start, uint32_t end) {
    auto n = end - start;
    bool swapped;
    for (uint32_t i = 0; i < n; i++) {
        swapped = false;
        for (uint32_t j = 0; j < n - i - 1; j++) {
            if (seq[start + j] > seq[start + j + 1]) {
                std::swap(seq[start + j], seq[start + j + 1]);
                swapped = true;
            }
        }
        if (!swapped) return;
    } 
}

// Simple Out-Of-Place Merge

inline void merge(parlay::slice<uint32_t*, uint32_t*> A, parlay::slice<uint32_t*, uint32_t*> B) {
    auto nA = A.size();
    auto nB = B.size();
    parlay::sequence<uint32_t> Aux(nA + nB);
    size_t i=0, j=0, k=0;
    while (i < nA && j < nB) {
        if (A[i] <= B[j]) {
        Aux[k++] = A[i++];
        } else {
        Aux[k++] = B[j++];
        }
    }

    while (i < nA) {
        Aux[k++] = A[i++];
    }
    while (j < nB) {
        Aux[k++] = B[j++];
    }

    for (size_t x = 0; x < nA; x++) {
        A[x] = Aux[x];
    }
    for (size_t x = 0; x < nB; x++) {
        B[x] = Aux[nA + x];
    }
}


inline void PairwiseSort(parlay::slice<uint32_t*, uint32_t*> block) {
  size_t n = block.size();
  // If n is odd, the last element has no partner
  size_t limit = n - (n % 2);

  // For each pair (2i, 2i+1), swap if out of order
  parlay::parallel_for(0, limit/2, [&](size_t i){
    size_t idx1 = 2*i;
    size_t idx2 = 2*i + 1;
    if (block[idx1] > block[idx2]) {
      std::swap(block[idx1], block[idx2]);
    }
  });
}

