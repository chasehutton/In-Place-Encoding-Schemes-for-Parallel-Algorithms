#pragma once

#include <cstdint>
#include <random>
#include <iostream>

#include "utils.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"

// Each SEGMENT can encode SEGMENT_SIZE / 2 bits
#define SEGMENT_SIZE 64

// Macros for Block j
// size b = constant * SEGMENT_SIZE
// Block j is encoded in S[j*b : j*b + b)

// T is encoded in S[j*b : j*b + SEGMENT_SIZE)
#define T_R(j) ReadBlock(seq, (j)*b, (j)*b + SEGMENT_SIZE)
#define T_W(j,v) WriteBlock(seq, (j)*b, (j)*b + SEGMENT_SIZE, v)

// inv is encoded in S[j*b + SEGMENT_SIZE : j*b + 2*SEGMENT_SIZE)
#define inv_R(j) ReadBlock(seq, (j)*b + SEGMENT_SIZE, (j)*b + 2*SEGMENT_SIZE)
#define inv_W(j,v) WriteBlock(seq, (j)*b + SEGMENT_SIZE, (j)*b + 2*SEGMENT_SIZE, v)

// R is encoded in S[j*b + 2*SEGMENT_SIZE : j*b + 3*SEGMENT_SIZE)
#define R_R(j) ReadBlock(seq, (j)*b + 2*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE)
#define R_W(j,v) WriteBlock(seq, (j)*b + 2*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE, v)

// C is encoded in S[j*b + 3*SEGMENT_SIZE :  j*b + 3*SEGMENT_SIZE + 2)
#define C_R(j) ReadBlock(seq, (j)*b + 3*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE + 2)
#define C_W(j,v) WriteBlock(seq, (j)*b + 3*SEGMENT_SIZE, (j)*b + 3*SEGMENT_SIZE + 2, v)

// D is encoded in S[j*b + 3*SEGMENT_SIZE + 2 : j*b + 3*SEGMENT_SIZE + 4)
#define D_R(j) ReadBlock(seq, (j)*b + 3*SEGMENT_SIZE + 2, (j)*b + 3*SEGMENT_SIZE + 4)
#define D_W(j,v) WriteBlock(seq, (j)*b + 3*SEGMENT_SIZE + 2, (j)*b + 3*SEGMENT_SIZE + 4, v)

// E is encoded in S[j*b + 3*SEGMENT_SIZE + 4 : j*b + 3*SEGMENT_SIZE + 6)
#define E_R(j) ReadBlock(seq, (j)*b + 3*SEGMENT_SIZE + 4, (j)*b + 3*SEGMENT_SIZE + 6)
#define E_W(j,v) WriteBlock(seq, (j)*b + 3*SEGMENT_SIZE + 4, (j)*b + 3*SEGMENT_SIZE + 6, v)

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
    auto nb = seq.size()/(2*b);

    parlay::parallel_for(0, nb, [&] (uint32_t i) {
        uint32_t low = 0;
        uint32_t high = nb;
        uint32_t mid;

        while (low < high) {
            mid = low + (high-low)/2;
            if (GET_ENDPOINT(i) > GET_ENDPOINT(mid + nb)) low = mid+1;
            else high = mid;
        }
        R_W(i, low);

        low = 0;
        high = nb;

        while (low < high) {
            mid = low + (high-low)/2;
            if (GET_ENDPOINT(i + nb) > GET_ENDPOINT(mid)) low = mid+1;
            else high = mid;
        }
        R_W(i + nb, low);
    });

    parlay::parallel_for(0, nb, [&] (uint32_t i) {
        T_W(i, R_R(i) + i);
        T_W(i + nb, R_R(i + nb) + i);
    });

    parlay::parallel_for(0, nb, [&] (uint32_t i) {
        uint32_t r1 = R_R(i);
        uint32_t r2 = R_R(i + nb);

        inv_W(i, T_R(nb + r1));
        inv_W(i + nb, T_R(r2));
    });
}

///////////////////////////////////////////////////////////////////////////////
// Determines when EndSort is finished. the size of seq is assumed to be     //
// divisible by and larger than b. b is assumed larger than SEGMENT_SIZE.      //
///////////////////////////////////////////////////////////////////////////////

inline bool Done(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {
    *flag = true;
    parlay::parallel_for(0, seq.size()/b, [&] (uint32_t i) {
        if (D_R(i) == 0) *flag = false;
    });

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
        });

        // 3) For i with E(i)=1 => do the swap
        parlay::parallel_for(0, nBlocks, [&] (uint32_t i) {
        if (D_R(i) == 0 && E_R(i) == 1) {
            D_W(i, 1);
            auto t = T_R(i);

            // Swap the blocks [i*b, i*b+b) and [t*b, t*b+b)
            SwapBlock(seq, i*b, i*b + b, t*b, t*b + b);

            // if T(i)==i => D(i)=1
            if (T_R(i) == i) {
                D_W(i, 1);
            }
        }
        });
    }
}

// Assumes nb >= 2
inline void Separate(parlay::sequence<uint32_t>& seq, uint32_t start, uint32_t end, uint32_t b) { 
    uint32_t nb = (end - start) / b;
    uint32_t i = (start / b) + (nb / 2) - 1;
    uint32_t j = std::min((start / b) + nb - 1, inv_R(i));

    // Save old inv
    uint32_t A_inv = inv_R(i);
    uint32_t B_inv = inv_R(j);

    // Make slices for blocks i and j
    // block i => [i*b, i*b + b)
    auto A = parlay::make_slice(seq.begin() + i*b,
                                seq.begin() + i*b + b);
    auto B = parlay::make_slice(seq.begin() + j*b,
                                seq.begin() + j*b + b);

    if (b <= 256) {
        BubbleSort(A, B);
    } else {
        PairwiseSort(A);
        PairwiseSort(B);
        merge(A, B);
    }

    // restore the inversion pointers
    inv_W(i, A_inv);
    inv_W(j, B_inv);
}

void SeqSort(parlay::sequence<uint32_t>& seq, uint32_t start, uint32_t end, uint32_t b) {
    assert(start % b == 0 && end % b == 0);
    uint32_t n = end - start;
    if (n <= 2*b) {
        if (b <= 512) {
            BubbleSort(seq, start, end);
        } else if (b <= 8192) {
            uint32_t half = start + (n/2);
            auto A = parlay::make_slice(seq.begin() + start,
                                        seq.begin() + half);
            auto B = parlay::make_slice(seq.begin() + half,
                                        seq.begin() + end);
                        
            PairwiseSort(A);
            PairwiseSort(B);
            merge(A, B);
        } else {
            // Parallel Merge Out-Of-Place
        }

        return;
    }

    uint32_t totalBlocks = n / b;   // # blocks in [start,end)
    uint32_t halfBlocks  = totalBlocks / 2;
    uint32_t mid = start + (halfBlocks * b);  // in multiples of b

    Separate(seq, start, end, b);

    parlay::par_do( 
        [&] {
            SeqSort(seq, start, mid, b);
        },
        [&] {
            SeqSort(seq, mid, end, b);
        }
    );
}

// GPT code to check if endsorted
bool IsEndSorted(const parlay::sequence<uint32_t>& seq, uint32_t b) {
  // number of blocks
  uint32_t nBlocks = seq.size() / b;
  
  // If only 0 or 1 block, trivially "end-sorted" by endpoints
  if (nBlocks <= 1) return true;

  // Check each adjacent pair of blocks
  for (uint32_t i = 1; i < nBlocks; i++) {
    // Compare block (i-1)'s endpoint with block i's endpoint
    if (GET_ENDPOINT(i - 1) > GET_ENDPOINT(i)) {
      std::cout << "Endpoints out of order between blocks "
                << (i - 1) << " and " << i
                << ": GET_ENDPOINT(" << (i - 1) << ")="
                << GET_ENDPOINT(i - 1)
                << " > GET_ENDPOINT(" << i << ")="
                << GET_ENDPOINT(i)
                << "\n";
      return false;
    }
  }

  return true;  // No violations found
}

// GPT code
bool CheckInversionPointers(parlay::sequence<uint32_t>& seq, uint32_t b) {
  // total #blocks across both halves
  uint32_t totalBlocks = seq.size() / b;
  // #blocks in each half
  uint32_t nb = totalBlocks / 2;  
  bool allGood = true;

  // ------------------------------------------------------------------------
  // 1) For blocks in the LEFT half [0..nb-1]:
  //    inv(i) should be in [nb.. 2*nb], 
  //    and be the first block in the right half whose endpoint >= GET_ENDPOINT(i).
  // ------------------------------------------------------------------------
  for (uint32_t i = 0; i < nb; i++) {
    // The endpoint of block i
    auto e_i = GET_ENDPOINT(i);

    // The "inversion pointer" for block i
    uint32_t inv_i = inv_R(i);

    // Check range. If we allow "inv_i = 2*nb" as a sentinel, that might be valid
    // depending on your code. Otherwise, typically inv_i < 2*nb.
    if (inv_i < nb || inv_i > 2*nb) {
      std::cerr << "ERROR: For block " << i << " in left half, inv_R(i)=" 
                << inv_i << " is out of [nb..2*nb)."
                << "\n";
      allGood = false;
      continue;
    }

    // 1a) For each block x in [nb.. inv_i),
    //     check GET_ENDPOINT(x) < e_i
    for (uint32_t x = nb; x < inv_i && x < 2*nb; x++) {
      if (GET_ENDPOINT(x) >= e_i) {
        std::cerr << "ERROR: block " << i << " endpoint=" << e_i
                  << " claims inv=" << inv_i
                  << " but block " << x
                  << " has endpoint=" << GET_ENDPOINT(x)
                  << " >= e_i.\n";
        allGood = false;
      }
    }

    // 1b) If inv_i < 2*nb, check GET_ENDPOINT(inv_i) >= e_i
    if (inv_i < 2*nb) {
      if (GET_ENDPOINT(inv_i) < e_i) {
        std::cerr << "ERROR: block " << i << " endpoint=" << e_i
                  << " claims inv=" << inv_i
                  << " but GET_ENDPOINT(inv_i)="
                  << GET_ENDPOINT(inv_i) << " < e_i.\n";
        allGood = false;
      }
    }
  }

  // ------------------------------------------------------------------------
  // 2) For blocks in the RIGHT half [nb..2*nb-1]:
  //    inv(i) should be in [0..nb], 
  //    and be the first block in the left half whose endpoint >= GET_ENDPOINT(i).
  // ------------------------------------------------------------------------
  for (uint32_t i = nb; i < 2*nb; i++) {
    auto e_i = GET_ENDPOINT(i);
    uint32_t inv_i = inv_R(i);

    // Check range: typically we expect inv_i < nb or inv_i=nb as a sentinel
    // if your code does that. Adjust if your logic differs (e.g. "inv_i<=nb").
    if (inv_i > nb) {
      std::cerr << "ERROR: For block " << i << " in right half, inv_R(i)=" 
                << inv_i << " is out of [0..nb].\n";
      allGood = false;
      continue;
    }

    // 2a) For each block x in [0.. inv_i),
    //     check GET_ENDPOINT(x) < e_i
    for (uint32_t x = 0; x < inv_i; x++) {
      if (GET_ENDPOINT(x) >= e_i) {
        std::cerr << "ERROR: block " << i << " endpoint=" << e_i
                  << " claims inv=" << inv_i
                  << " but block " << x
                  << " has endpoint=" << GET_ENDPOINT(x)
                  << " >= e_i.\n";
        allGood = false;
      }
    }

    // 2b) If inv_i < nb, check GET_ENDPOINT(inv_i) >= e_i
    if (inv_i < nb) {
      if (GET_ENDPOINT(inv_i) < e_i) {
        std::cerr << "ERROR: block " << i << " endpoint=" << e_i
                  << " claims inv=" << inv_i
                  << " but GET_ENDPOINT(inv_i)="
                  << GET_ENDPOINT(inv_i) << " < e_i.\n";
        allGood = false;
      }
    }
  }

  return allGood;
}

void Merge(parlay::sequence<uint32_t>& seq, uint32_t b) {
    assert(seq.size()/2 % b == 0);
    assert(seq.size() > b);
    assert(b % 2 == 0);
    assert(b >= 5*SEGMENT_SIZE);


    bool* flag = (bool*) std::malloc(sizeof(bool));
    *flag = false;

    std::cout << "Setting Up...\n\n";
    SetUp(seq, b);

    //assert(CheckInversionPointers(seq, b));

    std::cout << "End Sorting...\n\n";
    EndSort(seq, b, flag);

    assert(IsEndSorted(seq, b));

    std::cout << "Seq Sorting...\n\n";
    SeqSort(seq, 0, (uint32_t)seq.size(), b);

    std::free(flag);
}