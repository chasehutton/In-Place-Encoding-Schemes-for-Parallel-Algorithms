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
  if (b <= 4096) {
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
    std::array<uint32_t, 2*b> R;
    // parlay::sequence<uint32_t> R(2*b);
    parlay::internal::merge_into<parlay::copy_assign_tag>(A, B, parlay::make_slice(R),std::less<>());
    parlay::parallel_for(0, b, [&] (uint32_t i) {
        A[i] = R[i];
        B[i] = R[b + i];
    });
  }                
}

inline void reset_inv_encoding(parlay::sequence<uint32_t>& S, uint32_t start) {
    uint32_t i, idx1, idx2;

    for (i = 0; i < 32; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }
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

inline void reset_encoding(parlay::sequence<uint32_t>& S, uint32_t start) {
    uint32_t i, idx1, idx2;

    for (i = 0; i < 32; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start +  1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = 40; i < 72; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = 80; i < 81; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = bdiv2 + 9; i < bdiv2 + 41; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = bdiv2 + 49; i < bdiv2 + 113; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = bdiv2 + 121; i < 122; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }
}

inline void reset_encoding_minus_inv(parlay::sequence<uint32_t>& S, uint32_t start) {
    uint32_t i, idx1, idx2;

    for (i = 40; i < 72; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = 80; i < 81; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = bdiv2 + 9; i < bdiv2 + 41; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = bdiv2 + 49; i < bdiv2 + 113; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
    }

    for (i = bdiv2 + 121; i < 122; i++) {
        idx1 = 2*i + start;
        idx2 = 2*i + start + 1;
        if (S[idx1] > S[idx2]) {
            std::swap(S[idx1], S[idx2]);
        }
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

void ParReverseSplit(parlay::slice<uint32_t , uint32_t> A, parlay::slice<uint32_t , uint32_t> B) {
    auto Asize = A.size();
    auto Bsize = B.size();
    auto n = Asize + Bsize;
    if (n <= 2048) {
      for (int i = 0; i < std::floor((n) / 2); i++) {
        if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
        else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
        else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
      }
    } else {
      parlay::parallel_for(0, std::floor(n / 2), [&] (uint32_t i) {
        if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
        else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
        else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
      });
    }
  }
  
  void ParReverse(parlay::slice<uint32_t* , uint32_t*> A) {
    auto Asize = A.size();
    if (Asize <= 2048) {
      for (int i = 0; i < Asize / 2; i++) {
        std::swap(A[i], A[Asize - i - 1]);
      }
    } else {
      parlay::parallel_for(0, Asize / 2, [&] (uint32_t i) {
        std::swap(A[i], A[Asize - i - 1]);
      });
    }
  }

  void ParReverseSplit(parlay::slice<uint32_t *, uint32_t *> A, parlay::slice<uint32_t *, uint32_t *> B) {
    auto Asize = A.size();
    auto Bsize = B.size();
    auto n = Asize + Bsize;
    if (n <= 2048) {
      for (int i = 0; i < std::floor((n) / 2); i++) {
        if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
        else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
        else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
      }
    } else {
      parlay::parallel_for(0, std::floor(n / 2), [&] (uint32_t i) {
        if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
        else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
        else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
      });
    }
  }
  
  void SwapContiguousChunk(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t s, uint32_t m, uint32_t e) {
   int nA = A.size();
   int nB = B.size();
   // chunk contained in left side
   if (e < nA) {
       ParReverse(parlay::make_slice(A.begin()+s, A.begin()+e));
       ParReverse(parlay::make_slice(A.begin()+s, A.begin()+m));
       ParReverse(parlay::make_slice(A.begin()+m, A.begin()+e));
   }
   // chunk contained in right side
   else if (s >= nB) {
       ParReverse(parlay::make_slice(B.begin()+s-nA, B.begin()+e-nA));
       ParReverse(parlay::make_slice(B.begin()+s-nA, B.begin()+m-nA));
       ParReverse(parlay::make_slice(B.begin()+m-nA, B.begin()+e-nA));
   }
   // right side split
   else if (m < nA){
        ParReverseSplit(parlay::make_slice(A.begin()+s, A.end()), parlay::make_slice(B.begin(), B.begin()+e-nA));
        ParReverse(parlay::make_slice(A.begin()+s, A.begin()+m));
        ParReverseSplit(parlay::make_slice(A.begin()+m, A.end()), parlay::make_slice(B.begin(), B.begin()+e-nA));
    }
    // left side split
    else {
         ParReverseSplit(parlay::make_slice(A.begin()+s, A.end()), parlay::make_slice(B.begin(), B.begin()+e-nA));
         ParReverseSplit(parlay::make_slice(A.begin()+s, A.end()), parlay::make_slice(B.begin(), B.begin()+m-nA));
       ParReverse(parlay::make_slice(B.begin()+m-nA, B.begin()+e-nA));
    }
  }
