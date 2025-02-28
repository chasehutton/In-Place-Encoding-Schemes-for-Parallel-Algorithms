#pragma once

#include <cstdint>
#include <cstring> 

#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "parlay/primitives.h"
#include "block_size.h"

#define SEGMENT_SIZE 64
#define PADDING 16
static uint32_t SEGMENT_SIZE_t2 = 2*SEGMENT_SIZE; 
static uint32_t SEGMENT_SIZE_t3 = 3*SEGMENT_SIZE; 
static uint32_t SEGMENT_SIZE_t4 = 4*SEGMENT_SIZE; 
static uint32_t SEGMENT_SIZE_t5 = 5*SEGMENT_SIZE; 
static uint32_t PADDING_t2 = 2*PADDING;
static uint32_t PADDING_t3 = 3*PADDING;
static uint32_t PADDING_t4 = 4*PADDING;
static uint32_t PADDING_t5 = 5*PADDING;

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

uint32_t inline get_endpoint(parlay::sequence<uint32_t>& S, uint32_t j, int offset) {
    return S[j*b + b - 1 + offset];
}

inline void setup_end_sorted_position_buffer(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t v, int offset) {
   uint32_t seq_index = j*b + bdiv2 + PADDING + offset;
   write_block_64(S, seq_index + 2, S[seq_index]);
   S[seq_index] = v;
}

inline void restore_end_sorted_position_buffer(parlay::sequence<uint32_t>& S, uint32_t j, int offset) {
    uint32_t seq_index = j*b + bdiv2 + PADDING + offset;
    S[seq_index] = read_block_64(S, seq_index + 2);
}

uint32_t inline read_end_sorted_position(parlay::sequence<uint32_t>& S, uint32_t j, int offset) {
    uint32_t seq_index = b*j + bdiv2 + PADDING + offset;
    return read_block_64(S, seq_index);
}

void inline write_end_sorted_position(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, int offset) {
    uint32_t seq_index = b*j + bdiv2 + PADDING + offset;
    write_block_64(S, seq_index, value);
}

uint32_t inline read_inversion_pointer(parlay::sequence<uint32_t>& S, uint32_t j, int offset) {
    uint32_t seq_index = b*j + offset;
    return read_block_64(S, seq_index);
}

void inline write_inversion_pointer(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, int offset) {
    uint32_t seq_index = b*j + offset;
    write_block_64(S, seq_index, value);
}

uint32_t inline read_rank(parlay::sequence<uint32_t>& S, uint32_t j, int offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE + PADDING + offset;
    return read_block_64(S, seq_index);
}

void inline write_rank(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, int offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE + PADDING + offset;
    write_block_64(S, seq_index, value);
}

uint32_t inline read_coin_flip(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t k, int offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE + 2 + PADDING_t2 + 2*k + offset;
    return static_cast<uint32_t>(S[seq_index] > S[seq_index + 1]);
}

void inline write_coin_flips(parlay::sequence<uint32_t>& S, uint32_t j, uint64_t r, int offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE + 2 + PADDING_t2 + offset;
    write_block_64(S, seq_index, r);
}


uint32_t inline read_swap_flag(parlay::sequence<uint32_t>& S, uint32_t j, int offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE_t2 + PADDING_t2 + offset;
    return static_cast<uint32_t>(S[seq_index] > S[seq_index + 1]);
}

void inline write_swap_flag(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, int offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE_t2 + PADDING_t2 + offset;
    write_block_2(S, seq_index, value);
}


void inline mark_self(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t v, int offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE_t3 + 2 + PADDING_t3 + offset;
    write_block_2(S, seq_index, v);
}

uint32_t inline read_mark(parlay::sequence<uint32_t>& S, uint32_t j, int offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE_t3 + 2 + PADDING_t3 + offset;
    return static_cast<uint32_t>(S[seq_index] > S[seq_index + 1]); 
}


inline void swap_block_cpy_half(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B,
                         uint32_t block1_start, uint32_t block1_end,
                         uint32_t block2_start, uint32_t block2_end) {
    if (bdiv2 <= 4096) {
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
    
    if (b <= 4096) {
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

// phantom in S1
inline void swap_block_cpy(parlay::sequence<uint32_t>& S1, parlay::sequence<uint32_t>& S2, parlay::sequence<uint32_t>& phantom, 
                           uint32_t block1_start, uint32_t block1_end, uint32_t block2_start, uint32_t block2_end) {
    if (b <= 4096) {
        std::array<uint32_t, b> temp;

        // split block copied into tmp
        std::copy(phantom.begin(), phantom.end(), temp.begin());
        std::copy(S1.begin() + block1_start, S1.begin() + block1_end, temp.begin() + phantom.size());

        // S2 block copied into split block
        std::copy(S2.begin() + block2_start, S2.begin() + block2_start + phantom.size(), phantom.begin());
        std::copy(S2.begin() + block2_start + phantom.size(), S2.begin() + block2_end, S1.begin() + block1_start);

        // split block placed where S2 block was
        std::copy(temp.begin(), temp.end(), S2.begin() + block2_start);
    } else {

    }
}


inline void merge(parlay::slice<uint32_t*, uint32_t*> A,
                  parlay::slice<uint32_t*, uint32_t*> B) {
  uint32_t Asize = A.size();
  uint32_t Bsize = B.size();
  if (Asize == 0 || Bsize == 0) return;
  if (std::max(Asize, Bsize) <= 4096) {
    std::array<uint32_t, 2*b> temp;

    size_t i = 0; 
    size_t j = 0; 
    size_t k = 0;  

    while (i < Asize && j < Bsize) {
        if (A[i] <= B[j]) {
        temp[k++] = A[i++];
        } else {
        temp[k++] = B[j++];
        }
    }
    while (i < Asize) {
        temp[k++] = A[i++];
    }
    while (j < Bsize) {
        temp[k++] = B[j++];
    }

    for (size_t idx = 0; idx < Asize; idx++) {
        A[idx] = temp[idx];
    }
    for (size_t idx = 0; idx < Bsize; idx++) {
        B[idx] = temp[Asize + idx];
    }
  } else {
    std::array<uint32_t, 2*b> R;
    // parlay::sequence<uint32_t> R(2*b);
    parlay::internal::merge_into<parlay::copy_assign_tag>(A, B, parlay::make_slice(R),std::less<>());
    if (Asize <= 4096) {    
        for (size_t idx = 0; idx < Asize; idx++) {
            A[idx] = R[idx];
        }
        parlay::parallel_for(0, Bsize, [&] (uint32_t i) {
            B[i] = R[i + Bsize];
        });
    } else if (Bsize <= 4096) {
        for (size_t idx = 0; idx < Bsize; idx++) {
            B[idx] = R[idx + Asize];
        }
        parlay::parallel_for(0, Asize, [&] (uint32_t i) {
            A[i] = R[i];
        });

    } else {
        parlay::parallel_for(0, Asize, [&] (uint32_t i) {
            A[i] = R[i];
        });
        parlay::parallel_for(0, Bsize, [&] (uint32_t i) {
            B[i] = R[i + Bsize];
        });
    }   
  }                
}

inline void reset_inv_encoding(parlay::sequence<uint32_t>& S, uint32_t start, int offset) {
    write_inversion_pointer(S, start, 0, offset);
    write_rank(S, start, 0, offset);
    write_swap_flag(S, start, 0, offset);
    write_coin_flips(S, start, 0, offset);
    mark_self(S, start, 0, offset);
}

inline void reset_encoding_minus_inv(parlay::sequence<uint32_t>& S, uint32_t start, int offset) {
    write_rank(S, start, 0, offset);
    write_swap_flag(S, start, 0, offset);
    write_coin_flips(S, start, 0, offset);
    mark_self(S, start, 0, offset);
}

inline void reset_inv_encoding(parlay::slice<uint32_t*, uint32_t*> S) {
    uint32_t i, idx1, idx2;

    for (i = 0; i < 32; i++) {
        idx1 = 2*i;
        idx2 = 2*i + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }
}

// inline void reset_inv_encoding(parlay::sequence<uint32_t>& S, uint32_t start) {
//     write_inversion_pointer(S, start, 0);
//     write_rank(S, start, 0);
//     write_swap_flag(S, start, 0);
//     write_coin_flips(S, start, 0);
//     mark_self(S, start, 0);
// }

// inline void reset_encoding_minus_inv(parlay::sequence<uint32_t>& S, uint32_t start) {
//     write_rank(S, start, 0);
//     write_swap_flag(S, start, 0);
//     write_coin_flips(S, start, 0);
//     mark_self(S, start, 0);
// }


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

// void ParReverseSplit(parlay::slice<uint32_t , uint32_t> A, parlay::slice<uint32_t , uint32_t> B) {
//     auto Asize = A.size();
//     auto Bsize = B.size();
//     auto n = Asize + Bsize;
//     if (n <= 2048) {
//       for (int i = 0; i < std::floor((n) / 2); i++) {
//         if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
//         else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
//         else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
//       }
//     } else {
//       parlay::parallel_for(0, std::floor(n / 2), [&] (uint32_t i) {
//         if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
//         else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
//         else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
//       });
//     }
//   }
  
//   void ParReverse(parlay::slice<uint32_t* , uint32_t*> A) {
//     auto Asize = A.size();
//     if (Asize <= 2048) {
//       for (int i = 0; i < Asize / 2; i++) {
//         std::swap(A[i], A[Asize - i - 1]);
//       }
//     } else {
//       parlay::parallel_for(0, Asize / 2, [&] (uint32_t i) {
//         std::swap(A[i], A[Asize - i - 1]);
//       });
//     }
//   }

//   void ParReverseSplit(parlay::slice<uint32_t *, uint32_t *> A, parlay::slice<uint32_t *, uint32_t *> B) {
//     auto Asize = A.size();
//     auto Bsize = B.size();
//     auto n = Asize + Bsize;
//     if (n <= 2048) {
//       for (int i = 0; i < std::floor((n) / 2); i++) {
//         if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
//         else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
//         else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
//       }
//     } else {
//       parlay::parallel_for(0, std::floor(n / 2), [&] (uint32_t i) {
//         if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
//         else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
//         else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
//       });
//     }
//   }
  
//   void SwapContiguousChunk(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t s, uint32_t m, uint32_t e) {
//    int nA = A.size();
//    int nB = B.size();
//    // chunk contained in left side
//    if (e < nA) {
//        ParReverse(parlay::make_slice(A.begin()+s, A.begin()+e));
//        ParReverse(parlay::make_slice(A.begin()+s, A.begin()+m));
//        ParReverse(parlay::make_slice(A.begin()+m, A.begin()+e));
//    }
//    // chunk contained in right side
//    else if (s >= nB) {
//        ParReverse(parlay::make_slice(B.begin()+s-nA, B.begin()+e-nA));
//        ParReverse(parlay::make_slice(B.begin()+s-nA, B.begin()+m-nA));
//        ParReverse(parlay::make_slice(B.begin()+m-nA, B.begin()+e-nA));
//    }
//    // right side split
//    else if (m < nA){
//         ParReverseSplit(parlay::make_slice(A.begin()+s, A.end()), parlay::make_slice(B.begin(), B.begin()+e-nA));
//         ParReverse(parlay::make_slice(A.begin()+s, A.begin()+m));
//         ParReverseSplit(parlay::make_slice(A.begin()+m, A.end()), parlay::make_slice(B.begin(), B.begin()+e-nA));
//     }
//     // left side split
//     else {
//          ParReverseSplit(parlay::make_slice(A.begin()+s, A.end()), parlay::make_slice(B.begin(), B.begin()+e-nA));
//          ParReverseSplit(parlay::make_slice(A.begin()+s, A.end()), parlay::make_slice(B.begin(), B.begin()+m-nA));
//        ParReverse(parlay::make_slice(B.begin()+m-nA, B.begin()+e-nA));
//     }
//   }
