#pragma once

#include <cstdint>
#include <cstring> 

#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/utilities.h"
#include "parlay/primitives.h"

//thread_local static uint32_t swap_buffer[4096];
//thread_local static uint32_t merge_buffer[8192];


inline void SwapBlock(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
                      uint32_t block1_start, uint32_t block1_end,
                      uint32_t block2_start, uint32_t block2_end) {

    uint32_t size = block1_end - block1_start;
    for (int i = 0; i < size; i++) {
        std::swap(A[block1_start + i], B[block2_start + i]);
    }
}

inline void SwapBlock(parlay::slice<uint32_t*, uint32_t*> block1, 
                parlay::slice<uint32_t*, uint32_t*> block2) {
    
    uint32_t size = block1.size();
    for (int i = 0; i < size; i++) {
        std::swap(block1[i], block2[i]);
    }
}

inline void SwapBlockCpy(parlay::sequence<uint32_t>& seq, uint32_t block1_start, uint32_t block1_end,
                         uint32_t block2_start, uint32_t block2_end) {
    uint32_t size = block1_end - block1_start;

    if (size <= 2048) {
        parlay::sequence<uint32_t> temp(size);
        //auto temp = swap_buffer;
        std::copy(seq.begin() + block1_start, seq.begin() + block1_end, temp.begin());
        std::copy(seq.begin() + block2_start, seq.begin() + block2_end, seq.begin() + block1_start);
        std::copy(temp.begin(), temp.begin() + size, seq.begin() + block2_start);
  } else {
        parlay::parallel_for(0, size, [&] (uint32_t i) {
            std::swap(seq[block1_start + i], seq[block2_start + i]);
    });
  }
}

inline void SwapBlockCpy(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B,
                         uint32_t block1_start, uint32_t block1_end,
                         uint32_t block2_start, uint32_t block2_end) {
    uint32_t size = block1_end - block1_start;

    if (size <= 4096) {
        parlay::sequence<uint32_t> temp(size);
        std::copy(A.begin() + block1_start, A.begin() + block1_end, temp.begin());
        std::copy(B.begin() + block2_start, B.begin() + block2_end, A.begin() + block1_start);
        std::copy(temp.begin(), temp.begin() + size, B.begin() + block2_start);
  } else {
        parlay::parallel_for(0, size, [&] (uint32_t i) {
            std::swap(A[block1_start + i], B[block2_start + i]);
    });
  }
}

// inline void SwapBlockCpy(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B,
//                          uint32_t block1_start, uint32_t block1_end,
//                          uint32_t block2_start, uint32_t block2_end) {
//     uint32_t size = block1_end - block1_start;

//     if (size <= 4096) {
//         uint32_t temp[4096];
//         auto x = size*sizeof(uint32_t);
//         std::memmove(temp, A.data() + block1_start, x);
//         std::memmove(A.data() + block1_start, B.data() + block2_start, x);
//         std::memmove(B.data() + block2_start, temp, x);
//   } else {
//         parlay::parallel_for(0, size, [&] (uint32_t i) {
//             std::swap(A[block1_start + i], B[block2_start + i]);
//     });
//   }
// }


inline void SwapBlockCpy(parlay::slice<uint32_t*, uint32_t*> A, parlay::slice<uint32_t*, uint32_t*> B) {
    uint32_t size = A.size();   
    if (size <= 2048) {
        parlay::sequence<uint32_t> temp(size);
        std::copy(A.begin(), A.end(), temp.begin());
        std::copy(B.begin(), B.end(), A.begin());
        std::copy(temp.begin(), temp.begin() + size, B.begin());
    } else {
        parlay::parallel_for(0, size, [&] (uint32_t i) {
            std::swap(A[i], B[i]);
        });
    }
}

inline uint32_t ReadBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end) {
    //uint32_t size = end - start;
    uint32_t result = 0;
    uint8_t bit_pos = 0;
    for(int i = 0; i < 32; i++) {
        if (block[start + 2*i] > block[start + 2*i + 1]) result |= (1U << bit_pos);
        bit_pos++;
    }

    return result;
}

inline uint32_t ReadBlock(parlay::slice<uint32_t*, uint32_t*> block, uint32_t start, uint32_t end) {
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

inline void WriteBlock(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end, uint32_t value) {
    uint32_t size = end - start;
    for (int i = 0; i < size/2; i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = block[start + 2*i];
        uint32_t& second = block[start + 2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}

inline void WriteBlock128(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t end, uint64_t value) {
    uint32_t size = end - start;
    for (int i = 0; i < 64; i++) {
        bool bit = (value >> i) & 1U;
        uint32_t& first = block[start + 2*i];
        uint32_t& second = block[start + 2*i + 1];
        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        } 
    }
}

inline void WriteBlock2(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t value) {
    bool bit = value & 1U;
    uint32_t& first = block[start];
    uint32_t& second = block[start + 1];
    if ((!bit && first > second) || (bit && first < second)) {
        std::swap(first, second);
    }
}

inline void WriteBlock64(parlay::sequence<uint32_t>& block, uint32_t start, uint32_t value) {
    for (int i = 0; i < 32; ++i) {
        bool bit = (value >> i) & 1U;

        uint32_t& first = block[start + 2 * i];
        uint32_t& second = block[start + 2 * i + 1];

        if ((!bit && first > second) || (bit && first < second)) {
            std::swap(first, second);
        }
    }
}

inline void WriteBlock(parlay::slice<uint32_t*, uint32_t*> block, uint32_t start, uint32_t end, uint32_t value) {
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

// inline void merge(parlay::slice<uint32_t*, uint32_t*> A,
//                   parlay::slice<uint32_t*, uint32_t*> B) {
//   size_t n = A.size();
  
//   if (n <= 8192) {
//     uint32_t temp[8192];   

//     uint32_t* Aptr = A.begin();
//     uint32_t* Bptr = B.begin();
//     size_t i = 0, j = 0, k = 0;
    
//     while (i < n && j < n) {
//       if (Aptr[i] <= Bptr[j]) temp[k++] = Aptr[i++];
//       else                    temp[k++] = Bptr[j++];
//     }
//     while (i < n) temp[k++] = Aptr[i++];
//     while (j < n) temp[k++] = Bptr[j++];
//     auto x = n*sizeof(uint32_t);

//     std::memcpy(Aptr, temp, x);
//     std::memcpy(Bptr, temp + n, x);
//   }
//   else {
//     parlay::sequence<uint32_t> R(2 * n);
//     parlay::internal::merge_into<parlay::copy_assign_tag>(A, B, parlay::make_slice(R), std::less<>());
    
//     parlay::parallel_for(0, n, [&](size_t i) {
//       A[i] = R[i];
//       B[i] = R[n + i];
//     });
//   }
// }

inline void merge23(parlay::slice<uint32_t*, uint32_t*> A,
                  parlay::slice<uint32_t*, uint32_t*> B) {
  size_t n = A.size();


  if (n <= 8192) {
    parlay::sequence<uint32_t> temp(2 * n);

    size_t i = 0; 
    size_t j = 0; 
    size_t k = 0;  

    while (i < n && j < n) {
        if (A[i] <= B[j]) {
        temp[k++] = A[i++];
        } else {
        temp[k++] = B[j++];
        }
    }
    while (i < n) {
        temp[k++] = A[i++];
    }
    while (j < n) {
        temp[k++] = B[j++];
    }

    for (size_t idx = 0; idx < n; idx++) {
        A[idx] = temp[idx];
    }
    for (size_t idx = 0; idx < n; idx++) {
        B[idx] = temp[n + idx];
    }
  } else {
    parlay::sequence<uint32_t> R(n*2);
    parlay::internal::merge_into<parlay::copy_assign_tag>(A,B,parlay::make_slice(R),std::less<>());
    parlay::parallel_for(0, n, [&] (uint32_t i) {
        A[i] = R[i];
        B[i] = R[n + i];
    });
  }                
}

// // Simple Out-Of-Place Merge
// inline void merge(parlay::slice<uint32_t*, uint32_t*> A, parlay::slice<uint32_t*, uint32_t*> B) {
//     auto n = A.size();
//     parlay::sequence<uint32_t> temp(2 * n);

//     auto itA = A.begin();
//     auto itB = B.begin();
//     auto itTemp = temp.begin();

//     while (itA != A.end() && itB != B.end()) {
//         if (*itA <= *itB) {
//             *itTemp = *itA;
//             ++itA;
//         } else {
//             *itTemp = *itB;
//             ++itB;
//         }
//         ++itTemp;
//     }

//     std::copy(itA, A.end(), itTemp);
//     std::copy(itB, B.end(), itTemp);

//     std::copy(temp.begin(), temp.begin() + n, A.begin());
//     std::copy(temp.begin()+n, temp.end(), B.begin());
// }

inline void PairwiseSort(parlay::slice<uint32_t*, uint32_t*> block) {
  uint32_t n = block.size();
  uint32_t idx1 = 0;
  uint32_t idx2 = 0;

  if (n <= 8192) {
    for (int i = 0; i < n/2; i++) {
        idx1 = 2*i;
        idx2 = 2*i + 1;
        if (block[idx1] > block[idx2]) {
            std::swap(block[idx1], block[idx2]);
        }
    }
  } else {
    parlay::parallel_for(0, n/2, [&] (uint32_t i) {
        uint32_t idx1 = 2*i;
        uint32_t idx2 = 2*i + 1;
        if (block[idx1] > block[idx2]) {
            std::swap(block[idx1], block[idx2]);
        }
    });
  }
}

inline uint32_t ReadBlock64(parlay::sequence<uint32_t>& seq, uint32_t start) {
    uint32_t result = 0;

    for (int i = 0; i < 32; ++i) {
        result |= (static_cast<uint32_t>(seq[start + 2 * i] > seq[start + 2 * i + 1]) << i);
    }

    return result;
}

