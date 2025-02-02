#pragma once

#include <cstdint>
#include <random>
#include <iostream>
#include <chrono>


#include "utils.h"
#include "parlay/sequence.h"
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/utilities.h"
#include "parlay/internal/binary_search.h"

#define OPTIMIZE 0

// Each SEGMENT can encode SEGMENT_SIZE / 2 bits
#define SEGMENT_SIZE 64

// Macros for Block j
// size b = constant * SEGMENT_SIZE
// Block j is encoded in S[j*b : j*b + b)

// T is encoded in S[j*b : j*b + SEGMENT_SIZE)
// #define T_R(j) ReadBlock64(seq, (j)*b)
#define T_R(j) ReadBlock(seq, (j)*b, (j)*b + SEGMENT_SIZE)
#define T_W(j,v) WriteBlock(seq, (j)*b, (j)*b + SEGMENT_SIZE, v)

// inv is encoded in S[j*b + SEGMENT_SIZE : j*b + 2*SEGMENT_SIZE)
// #define inv_R(j) ReadBlock64(seq, (j)*b + SEGMENT_SIZE)
#define inv_R(j) ReadBlock(seq, (j)*b + SEGMENT_SIZE, (j)*b + 2*SEGMENT_SIZE)
#define inv_W(j,v) WriteBlock(seq, (j)*b + SEGMENT_SIZE, (j)*b + 2*SEGMENT_SIZE, v)

// R is encoded in S[j*b + 2*SEGMENT_SIZE : j*b + 3*SEGMENT_SIZE)
// #define R_R(j) ReadBlock64(seq, (j)*b + 2*SEGMENT_SIZE)
#define R_R(j) ReadBlock(seq, (j)*b + 2*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE)
#define R_W(j,v) WriteBlock(seq, (j)*b + 2*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE, v)

// C is encoded in S[j*b + 3*SEGMENT_SIZE :  j*b + 3*SEGMENT_SIZE + 2)
#define C_R(j) (static_cast<uint32_t>(seq[(j)*b + 3*SEGMENT_SIZE] > seq[(j)*b + 3*SEGMENT_SIZE + 1]))
// #define C_R(j) ReadBlock(seq, (j)*b + 3*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE + 2)
#define C_W(j,v) WriteBlock(seq, (j)*b + 3*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE + 2, v)

// D is encoded in S[j*b + 3*SEGMENT_SIZE + 2 : j*b + 3*SEGMENT_SIZE + 4)
#define D_R(j) (static_cast<uint32_t>(seq[(j)*b + 3*SEGMENT_SIZE + 2] > seq[(j)*b + 3*SEGMENT_SIZE + 3]))
//#define D_R(j) ReadBlock(seq, (j)*b + 3*SEGMENT_SIZE + 2, (j)*b + 3*SEGMENT_SIZE + 4)
#define D_W(j,v) WriteBlock(seq, (j)*b + 3*SEGMENT_SIZE + 2, (j)*b + 3*SEGMENT_SIZE + 4, v)

// E is encoded in S[j*b + 3*SEGMENT_SIZE + 4 : j*b + 3*SEGMENT_SIZE + 6)
#define E_R(j) (static_cast<uint32_t>(seq[(j)*b + 3*SEGMENT_SIZE + 4] > seq[(j)*b + 3*SEGMENT_SIZE + 5]))
// #define E_R(j) ReadBlock(seq, (j)*b + 3*SEGMENT_SIZE + 4, (j)*b + 3*SEGMENT_SIZE + 6)
#define E_W(j,v) WriteBlock(seq, (j)*b + 3*SEGMENT_SIZE + 4, (j)*b + 3*SEGMENT_SIZE + 6, v)

#define T2_R(j) ReadBlock(seq, (j)*b + (b/2), (j)*b + (b/2) +  SEGMENT_SIZE)
#define T2_W(j,v) WriteBlock(seq, (j)*b + (b/2), (j)*b + (b/2) +  SEGMENT_SIZE, v)

#define inv2_R(j) ReadBlock(seq, (j)*b + (b/2) + SEGMENT_SIZE, (j)*b + (b/2) + 2*SEGMENT_SIZE)
#define inv2_W(j,v) WriteBlock(seq, (j)*b + (b/2) + SEGMENT_SIZE, (j)*b + (b/2) + 2*SEGMENT_SIZE, v)

#define R2_R(j) ReadBlock(seq, (j)*b + (b/2) + 2*SEGMENT_SIZE, (j)*b + (b/2) + 3*SEGMENT_SIZE)
#define R2_W(j,v) WriteBlock(seq, (j)*b + (b/2) + 2*SEGMENT_SIZE, (j)*b + (b/2) + 3*SEGMENT_SIZE, v)

#define C2_R(j) (static_cast<uint32_t>(seq[(j)*b + (b/2) + 3*SEGMENT_SIZE] > seq[(j)*b + (b/2) + 3*SEGMENT_SIZE + 1]))
#define C2_W(j,v) WriteBlock(seq, (j)*b + (b/2) + 3*SEGMENT_SIZE, (j)*b + (b/2) + 3*SEGMENT_SIZE + 2, v)

#define D2_R(j) (static_cast<uint32_t>(seq[(j)*b + (b/2) + 3*SEGMENT_SIZE + 2] > seq[(j)*b + (b/2) + 3*SEGMENT_SIZE + 3]))
#define D2_W(j,v) WriteBlock(seq, (j)*b + (b/2) + 2, (j)*b + (b/2) + 4, v)

#define E2_R(j) (static_cast<uint32_t>(seq[(j)*b + (b/2) + 3*SEGMENT_SIZE + 4] > seq[(j)*b + 3*SEGMENT_SIZE + (b/2) + 5]))
#define E2_W(j,v) WriteBlock(seq, (j)*b + (b/2) + 3*SEGMENT_SIZE + 4, (j)*b + (b/2) + 3*SEGMENT_SIZE + 6, v)



// The "endpoint" of block j is the element at index (j*b + b -1).
// That is the last element in the half-open block S[j*b : j*b + b).
#define GET_ENDPOINT(j) seq[(j)*b + (b - 1)]

///////////////////////////////////////////////////////////////////////////////////////////
// Each worker thread gets a workspace of size c*b where b is the given non-encoding     //
// block size and c is a constant to be determined (currently set to be 4).              //
// The first 3*b uint32_t values of each workspace are reserved for the out-of-place     //
// sequential merging done in the Separate function. The remaining values are used       //
// for any miscellaneous auxilary storage.                                               //
///////////////////////////////////////////////////////////////////////////////////////////

static std::vector<uint32_t*> workspaces;

void InitWorkspaces(std::size_t n) {
    auto nw = parlay::num_workers();
    workspaces.resize(nw);
    
    for (int i = 0; i < nw; i++) {
        workspaces[i] = (uint32_t*) std::malloc(n * sizeof(uint32_t));
    }
}

void FreeWorkspaces() {
    for (auto* ptr : workspaces) {
        std::free(ptr);
    }

    workspaces.clear();
}

///////////////////////////////////////////////////////////////////////////
// Computes and encodes the end-sorted position and inversion pointers   //
// for each block.                                                       //
///////////////////////////////////////////////////////////////////////////

void SetUp(parlay::sequence<uint32_t>& seq, uint32_t b) {
    // # of blocks in each half
    auto halfBlocks = seq.size()/(2*b);

    parlay::parallel_for(0, halfBlocks, [&] (uint32_t i) {
        uint32_t low = 0;
        uint32_t high = halfBlocks;
        uint32_t mid;

        while (low < high) {
            mid = low + (high-low)/2;
            if (GET_ENDPOINT(i) > GET_ENDPOINT(mid + halfBlocks)) low = mid+1;
            else high = mid;
        }
        // if (low == halfBlocks + 1) low = halfBlocks;
        R_W(i, low);

        low = 0;
        high = halfBlocks;

        while (low < high) {
            mid = low + (high-low)/2;
            if (GET_ENDPOINT(i + halfBlocks) > GET_ENDPOINT(mid)) low = mid+1;
            else high = mid;
        }
        // if (low == halfBlocks + 1) low = halfBlocks;
        R_W(i + halfBlocks, low);
    });

    parlay::parallel_for(0, halfBlocks, [&] (uint32_t i) {
        T_W(i, R_R(i) + i);
        T_W(i + halfBlocks, R_R(i + halfBlocks) + i);
    });

    parlay::parallel_for(0, halfBlocks, [&] (uint32_t i) {
        uint32_t r1 = R_R(i);
        uint32_t r2 = R_R(i + halfBlocks);

        uint32_t v1 = (r1 >= halfBlocks) ? seq.size() : T_R(halfBlocks + r1);
        uint32_t v2 = (r2 >= halfBlocks) ? seq.size() : T_R(r2);

        inv_W(i, v1);
        inv_W(i + halfBlocks, v2);
    });
}

///////////////////////////////////////////////////////////////////////////////
// Determines when EndSort is finished. the size of seq is assumed to be     //
// divisible by and larger than b. b is assumed larger than SEGMENT_SIZE.      //
///////////////////////////////////////////////////////////////////////////////

inline bool Done(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {
    std::atomic<bool> done(true);
    auto nBlocks = seq.size() / b;
    parlay::parallel_for(0, nBlocks, [&](uint32_t i) {
        if (D_R(i) == 0) done.store(false, std::memory_order_relaxed);
    });
    *flag = done.load(std::memory_order_relaxed);
    
    //*flag = true;
    // parlay::parallel_for(0, seq.size()/b, [&] (uint32_t i) {
    //     if (D_R(i) == 0) *flag = false;
    // });
    return *flag;
}

void EndSort(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {
    uint32_t nBlocks = seq.size() / b;

    // Initialize E=0, D=1 if T(i)=i else 0
    parlay::parallel_for(0, nBlocks, [&] (uint32_t i) {
        E_W(i, 0);
        if (T_R(i) == i) {
            D_W(i, 1);
        } else {
            D_W(i, 0);
        }
    });

    int itcount = 0;
    while (!Done(seq, b, flag)) {
        itcount++;
        // if (itcount > 100000) assert(false && "too long");

        // 1) Flip a random coin for each block i that is undone => if D(i)=0
        parlay::parallel_for(0, nBlocks, [&] (uint32_t i) {
        if (D_R(i) == 0) {
            auto r = parlay::random(131542391u 
                                    + itcount * 0x9e3779b97f4a7c15ULL 
                                    + i);
            uint8_t bit = static_cast<uint8_t>(r.ith_rand(0) & 1ULL);
            C_W(i, bit);
        }
        });

        // 2) If i is undone, C(i)=1, and C(T(i))=0 => E(i)=1
        parlay::parallel_for(0, nBlocks, [&] (uint32_t i) {
            if (D_R(i) == 0 && C_R(i) == 1 && C_R( T_R(i) ) == 0) {
                E_W(i, 1);
            }

            T2_W(i, T_R(i));
            E2_W(i, E_R(i));
        });

        parlay::parallel_for(0, nBlocks, [&] (uint32_t i) {
            if (E2_R(i) == 1) {
                D_W(i, 1);
                R_W(i, i);
                auto t = T2_R(i);

                // Swap half blocks
                SwapBlockCpy(seq, i*b, i*b + b/2, t*b, t*b + b/2);
            }
        });

        parlay::parallel_for(0, nBlocks, [&] (uint32_t i) {
            if (E_R(i) == 1) {
                E_W(i, 0);
                auto t = R_R(i);
                if (T_R(t) == t) {
                    D_W(t, 1);
                }
                SwapBlockCpy(seq, i*b + b/2, i*b + b, t*b + b/2, t*b + b);
            }
        });

        // for (int i = 0; i < nBlocks; i++) {
        //     if (T_R(i) == i) {
        //         D_W(i, 1);
        //     }
        // }
        // for(int i = 0; i < nBlocks; i++) {
        // if (D_R(i) == 0 && E_R(i) == 1) {
        //     D_W(i, 1);
        //     auto t = T_R(i);

        //     // Swap the blocks [i*b, i*b+b) and [t*b, t*b+b)
        //     SwapBlock(seq, i*b, i*b + b, t*b, t*b + b);

        //     // if T(i)==i => D(i)=1
        //     if (T_R(i) == i) {
        //         D_W(i, 1);
        //     }
        // }
        // };
    }
}

bool CheckSorted(parlay::sequence<uint32_t>& seq) {
    bool x = true;
    for (int i = 0; i < seq.size() - 1; i++) {
        if (seq[i] > seq[i+1]) {
            x = false;
        }
    }

    return x;
}


bool CheckSorted(parlay::slice<uint32_t*, uint32_t*> A) {
    bool x = true;
    for (int i = 0; i < A.size() - 1; i++) {
        if (A[i] > A[i+1]) {
            x = false;
        }
    }

    return x;
}

// Assumes nb >= 2
inline void Separate(parlay::sequence<uint32_t>& seq, uint32_t start, uint32_t end, uint32_t b, bool base_case) { 
    uint32_t nBlocks = (end - start) / b;
    uint32_t i = (start / b) + (nBlocks / 2) - 1;
    // uint32_t j = std::min(inv_R(i), start/b + nBlocks - 1);

    // assert(start % b == 0);
    // assert(end % b == 0);
    // assert((end - start) % b == 0);

    uint32_t e = start/b + nBlocks - 1;
    uint32_t i_inv = inv_R(i);

    if (i_inv > e && !base_case) return;

    uint32_t j = std::min(i_inv, e);

    // Save old inv
    uint32_t A_inv = inv_R(i);
    uint32_t B_inv = inv_R(j);

    // Make slices for blocks i and j
    // block i => [i*b, i*b + b)
    auto A = parlay::make_slice(seq.begin() + i*b,
                                seq.begin() + i*b + b);
    auto B = parlay::make_slice(seq.begin() + j*b,
                                seq.begin() + j*b + b);

    if (b <= 128) {
       BubbleSort(A, B);  
    //    parlay::integer_sort_inplace(parlay::make_slice(seq.begin() + start, seq.begin() + end)); 
    } else {
        PairwiseSort(A);
        PairwiseSort(B);
        // if (!CheckSorted(A)) {
        //     for (int k = 0; k < A.size(); k++) {
        //         std::cout << "\n" << A[k] << "\n";
        //         std::cout.flush();
        //     }
        //     assert(false);
        // } else if (!CheckSorted(B)) {
        //     for (int k = 0; k < B.size(); k++) {
        //         std::cout << "\n" << B[k] << "\n";
        //         std::cout.flush();
        //     }
        //     assert(false);
        // }
        // assert(CheckSorted(A));
        // assert(CheckSorted(B));
        merge(A, B);
    }

    if (!base_case) {
        inv_W(i, A_inv);
        inv_W(j, B_inv);
    }
}

void SeqSort(parlay::sequence<uint32_t>& seq, uint32_t start, uint32_t end, uint32_t b) {
    //assert(start % b == 0 && end % b == 0);
    uint32_t n = end - start;
    uint32_t nBlocks = n / b;   // # blocks in [start,end)
    uint32_t halfBlocks  = nBlocks / 2;
    uint32_t mid = start + (halfBlocks * b);  // in multiples of b

    if (OPTIMIZE) {

    } else { 
        if (n == 2*b) {
            Separate(seq, start, end, b, true);
            return;
        } else if (n == b) {
            PairwiseSort(parlay::make_slice(seq.begin() + start, seq.begin() + end));
            return;
        }
    }
    Separate(seq, start, end, b, false);

    parlay::par_do( 
        [&] {
            SeqSort(seq, start, mid, b);
        },
        [&] {
            SeqSort(seq, mid, end, b);
        }
    );
}

bool CheckInversionPointers(parlay::sequence<uint32_t>& seq, uint32_t b) {
    uint32_t nb = seq.size()/(2*b);
    bool x = true;
    for (int i = 0; i < nb; i++){
        uint32_t e_i = GET_ENDPOINT(i);
        uint32_t inv_i = R_R(i);

        for (int j = nb; j < 2*nb; j++) {
            uint32_t t = GET_ENDPOINT(j);
            if (j < inv_i + nb) {
                if (e_i < t) {
                    std::cout << "Block " << i << " claims inversion with block " << inv_i + nb << " but has inversion with block " << j << "\n";
                    x = false;
                }
            } else if (j > inv_i + nb) {
                if (e_i > t) {
                    std::cout << "Block " << i << " claims inversion with block " << inv_i + nb << " but has inversion with block " << j << "\n";
                    x = false;
                }
            }
        }

        e_i = GET_ENDPOINT(i + nb);
        inv_i = R_R(i + nb);

        for (int j = 0; j < nb; j++) {
            uint32_t t = GET_ENDPOINT(j);
            if (j < inv_i) {
                if (e_i < t) {
                    std::cout << "Block " << i + nb << " claims inversion with block " << inv_i << " but has inversion with block " << j << "\n";
                    x = false;
                }
            } else if (j > inv_i) {
                if (e_i > t) {
                    std::cout << "Block " << i + nb << " claims inversion with block " << inv_i << " but has inversion with block " << j << "\n";
                    x = false;
                }
            }
        }
    }

    return x;
}

bool CheckEndSorted(parlay::sequence<uint32_t>& seq, uint32_t b) {
    bool x = true;
    for (int i = 0; i < seq.size() / b - 1; i++) {
        if (GET_ENDPOINT(i) > GET_ENDPOINT(i+1)) {
            x = false;
        }
    }

    return x;
}



void Merge(parlay::sequence<uint32_t>& seq, uint32_t b) {
    assert(seq.size()/2 % b == 0);
    assert(seq.size()/b % 2 == 0);
    assert(seq.size() > b);
    assert(b % 2 == 0);
    assert(b >= 5*SEGMENT_SIZE);

    bool* flag = (bool*) std::malloc(sizeof(bool));
    *flag = false;
    // std::cout << "Setting Up...\n\n";
    auto start = std::chrono::high_resolution_clock().now();
    SetUp(seq, b);
    auto end = std::chrono::high_resolution_clock().now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    //assert(CheckInversionPointers(seq, b));
    std::cout << "Time for SetUp: " << time.count() << " \n";
    // std::cout << "End Sorting...\n\n";
    start = std::chrono::high_resolution_clock().now();
    EndSort(seq, b, flag);
    end = std::chrono::high_resolution_clock().now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    // //assert(CheckEndSorted(seq, b));
    std::cout << "Time for EndSort: " << time.count() << " \n";
    // std::cout << "Seq Sorting...\n\n";
    start = std::chrono::high_resolution_clock().now();
    SeqSort(seq, 0, (uint32_t)seq.size(), b);
    end = std::chrono::high_resolution_clock().now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "Time for SeqSort: " << time.count() << " \n";

    // for (int i = 0; i < seq.size()/b; i++) {
    //     parlay::integer_sort_inplace(parlay::make_slice(seq.begin() + i*b, seq.begin() + (i+1)*b));
    // }

    // if (!CheckSorted(seq)) {
    //     std::cout << "\n\nNot Sorted. Value of rightmost left half block endpoint: " << t << "\n\n"; 
    // }
    //assert(CheckSorted(seq));

    std::free(flag);
}
