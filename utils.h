#pragma once

#include <cstdint>
#include <cstring> 

#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "parlay/primitives.h"
#include "block_size.h"

inline void swap_block_cpy_half(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B,
                         uint32_t block1_start, uint32_t block1_end,
                         uint32_t block2_start, uint32_t block2_end) {
    if (bdiv2 <= 8192) {
        std::array<uint32_t, bdiv2> temp;
        //parlay::sequence<uint32_t> temp(size);
        std::copy(A.begin() + block1_start, A.begin() + block1_end, temp.begin());
        std::copy(B.begin() + block2_start, B.begin() + block2_end, A.begin() + block1_start);
        std::copy(temp.begin(), temp.begin() + bdiv2, B.begin() + block2_start);
  } else {
        parlay::parallel_for(0, bdiv2, [&] (uint32_t i) {
            std::swap(A[block1_start + i], B[block2_start + i]);
    });
  }
}

inline void swap_block_cpy(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B,
    uint32_t block1_start, uint32_t block1_end,
    uint32_t block2_start, uint32_t block2_end) {
    
    if (b <= 8192) {
        std::array<uint32_t, b> temp;
        std::copy(A.begin() + block1_start, A.begin() + block1_end, temp.begin());
        std::copy(B.begin() + block2_start, B.begin() + block2_end, A.begin() + block1_start);
        std::copy(temp.begin(), temp.begin() + b, B.begin() + block2_start);
    } else {
        parlay::parallel_for(0, b, [&] (uint32_t i) {
        std::swap(A[block1_start + i], B[block2_start + i]);
    });
    }
}

inline uint32_t read_block_64(parlay::sequence<uint32_t>& S, uint32_t start) {
    uint32_t result = 0;
    uint32_t bit_pos = 0;
    for (int i = 0; i < 32; i++) {
        if (S[start + 2*i] > S[start + 2*i + 1]) result |= (1U << bit_pos);
        bit_pos++;
    }

    return result;
}

inline uint32_t read_block_128(parlay::sequence<uint32_t>& S, uint32_t start) {
    uint32_t result = 0;
    uint32_t bit_pos = 0;
    for (int i = 0; i < 64; i++) {
        if (S[start + 2*i] > S[start + 2*i + 1]) result |= (1U << bit_pos);
        bit_pos++;
    }

    return result;
}

inline void write_block_2(parlay::sequence<uint32_t>& S, uint32_t start, uint32_t value) {
    bool bit = value & 1U;
    uint32_t& first = S[start];
    uint32_t& second = S[start + 1];
    if ((!bit && first > second) || (bit && first < second)) {
        std::swap(first, second);
    }
}

inline void write_block_64(parlay::sequence<uint32_t>& S, uint32_t start, uint64_t value) {
    for (int i = 0; i < 32; i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = S[start + 2*i];
        uint32_t& second = S[start + 2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}

inline void write_block_128(parlay::sequence<uint32_t>& S, uint32_t start, uint64_t value) {
    for (int i = 0; i < 64; i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = S[start + 2*i];
        uint32_t& second = S[start + 2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}

inline void merge(parlay::slice<uint32_t*, uint32_t*> A,
                  parlay::slice<uint32_t*, uint32_t*> B) {
  if (b <= 8192) {
    std::array<uint32_t, 2*b> temp;

    size_t i = 0; 
    size_t j = 0; 
    size_t k = 0;  

    while (i < b && j < b) {
        if (A[i] <= B[j]) {
        temp[k++] = A[i++];
        } else {
        temp[k++] = B[j++];
        }
    }
    while (i < b) {
        temp[k++] = A[i++];
    }
    while (j < b) {
        temp[k++] = B[j++];
    }

    for (size_t idx = 0; idx < b; idx++) {
        A[idx] = temp[idx];
    }
    for (size_t idx = 0; idx < b; idx++) {
        B[idx] = temp[b + idx];
    }
  } else {
    // std::array<uint32_t, 2*b> R;
    parlay::sequence<uint32_t> R(2*b);
    parlay::internal::merge_into<parlay::copy_assign_tag>(A, B, parlay::make_slice(R),std::less<>());
    parlay::parallel_for(0, b, [&] (uint32_t i) {
        A[i] = R[i];
        B[i] = R[b + i];
    });
  }                
}

inline void pairwise_sort(parlay::slice<uint32_t*, uint32_t*> S) {
  uint32_t idx1 = 0;
  uint32_t idx2 = 0;

  if (b <= 8192) {
    for (int i = 0; i < bdiv2; i++) {
        idx1 = 2*i;
        idx2 = 2*i + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }
  } else {
    parlay::parallel_for(0, bdiv2, [&] (uint32_t i) {
        uint32_t idx1 = 2*i;
        uint32_t idx2 = 2*i + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    });
  }
}
