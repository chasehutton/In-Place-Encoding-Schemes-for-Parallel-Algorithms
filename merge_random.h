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

#define SEGMENT_SIZE 64
static uint32_t SEGMENT_SIZE_t2 = 2*SEGMENT_SIZE; 
static uint32_t SEGMENT_SIZE_t3 = 3*SEGMENT_SIZE; 
static uint32_t SEGMENT_SIZE_t4 = 4*SEGMENT_SIZE; 
static uint32_t SEGMENT_SIZE_t5 = 5*SEGMENT_SIZE; 
uint32_t bdiv2 = 0;

#define TA_R(j) ReadBlock(A, (j)*b, (j)*b + SEGMENT_SIZE)
#define TA_W(j,v) WriteBlock(A, (j)*b, (j)*b + SEGMENT_SIZE, v)

#define invA_R(j) ReadBlock(A, (j)*b + SEGMENT_SIZE, (j)*b + SEGMENT_SIZE_t2)
#define invA_W(j,v) WriteBlock(A, (j)*b + SEGMENT_SIZE, (j)*b + SEGMENT_SIZE_t2, v)

#define RA_R(j) ReadBlock(A, (j)*b + SEGMENT_SIZE_t2, (j)*b + SEGMENT_SIZE_t3)
#define RA_W(j,v) WriteBlock(A, (j)*b + SEGMENT_SIZE_t2, (j)*b + SEGMENT_SIZE_t3, v)

#define CA_R(j, k) (static_cast<uint32_t>(A[(j)*b + SEGMENT_SIZE_t3 + 2*(k)] > A[(j)*b + SEGMENT_SIZE_t3 + 2*(k) + 1]))
#define CA_W(j,v) WriteBlock128(A, (j)*b + SEGMENT_SIZE_t3, (j)*b + SEGMENT_SIZE_t5, v)

#define DA_R(j) (static_cast<uint32_t>(A[(j)*b + SEGMENT_SIZE_t5] > A[(j)*b + SEGMENT_SIZE_t5 + 1]))
#define DA_W(j,v) WriteBlock(A, (j)*b + SEGMENT_SIZE_t5, (j)*b + SEGMENT_SIZE_t5 + 2, v)

#define EA_R(j) (static_cast<uint32_t>(A[(j)*b + SEGMENT_SIZE_t5 + 2] > A[(j)*b + SEGMENT_SIZE_t5 + 3]))
#define EA_W(j,v) WriteBlock(A, (j)*b + SEGMENT_SIZE_t5 + 2, (j)*b + SEGMENT_SIZE_t5 + 4, v)

#define CAS_R(j) (static_cast<uint32_t>(A[(j)*b + SEGMENT_SIZE_t5 + 4] > A[(j)*b + SEGMENT_SIZE_t5 + 5]))
#define CAS_W(j,v) WriteBlock(A, (j)*b + SEGMENT_SIZE_t5 + 4, (j)*b + SEGMENT_SIZE_t5 + 6, v)

#define TA2_R(j) ReadBlock(A, (j)*b + bdiv2, (j)*b + bdiv2+  SEGMENT_SIZE)
#define TA2_W(j,v) WriteBlock(A, (j)*b + bdiv2, (j)*b + bdiv2 +  SEGMENT_SIZE, v)

#define EA2_R(j) (static_cast<uint32_t>(A[(j)*b + (bdiv2) + SEGMENT_SIZE_t3 + 4] > A[(j)*b + SEGMENT_SIZE_t3 + (bdiv2) + 5]))
#define EA2_W(j,v) WriteBlock(A, (j)*b + (bdiv2) + SEGMENT_SIZE_t3 + 4, (j)*b + (bdiv2) + SEGMENT_SIZE_t3 + 6, v)

#define TB_R(j) ReadBlock(B, (j)*b, (j)*b + SEGMENT_SIZE)
#define TB_W(j,v) WriteBlock(B, (j)*b, (j)*b + SEGMENT_SIZE, v)

#define invB_R(j) ReadBlock(B, (j)*b + SEGMENT_SIZE, (j)*b + SEGMENT_SIZE_t2)
#define invB_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE, (j)*b + SEGMENT_SIZE_t2, v)

#define RB_R(j) ReadBlock(B, (j)*b + SEGMENT_SIZE_t2, (j)*b + SEGMENT_SIZE_t3)
#define RB_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE_t2, (j)*b + SEGMENT_SIZE_t3, v)

#define CB_R(j, k) (static_cast<uint32_t>(B[(j)*b + SEGMENT_SIZE_t3 + 2*(k)] > B[(j)*b + SEGMENT_SIZE_t3 + 2*(k) + 1]))
#define CB_W(j,v) WriteBlock128(B, (j)*b + SEGMENT_SIZE_t3, (j)*b + SEGMENT_SIZE_t5, v)

// #define CB_R(j) (static_cast<uint32_t>(B[(j)*b + SEGMENT_SIZE_t3] > B[(j)*b + SEGMENT_SIZE_t3 + 1]))
// #define CB_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE_t3, (j)*b + SEGMENT_SIZE_t3 + 2, v)

// #define DB_R(j) (static_cast<uint32_t>(B[(j)*b + SEGMENT_SIZE_t3 + 2] > B[(j)*b + SEGMENT_SIZE_t3 + 3]))
// #define DB_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE_t3 + 2, (j)*b + SEGMENT_SIZE_t3 + 4, v)

// #define EB_R(j) (static_cast<uint32_t>(B[(j)*b + SEGMENT_SIZE_t3 + 4] > B[(j)*b + SEGMENT_SIZE_t3 + 5]))
// #define EB_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE_t3 + 4, (j)*b + SEGMENT_SIZE_t3 + 6, v)

#define DB_R(j) (static_cast<uint32_t>(B[(j)*b + SEGMENT_SIZE_t5] > B[(j)*b + SEGMENT_SIZE_t5 + 1]))
#define DB_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE_t5, (j)*b + SEGMENT_SIZE_t5 + 2, v)

#define EB_R(j) (static_cast<uint32_t>(B[(j)*b + SEGMENT_SIZE_t5 + 2] > B[(j)*b + SEGMENT_SIZE_t5 + 3]))
#define EB_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE_t5 + 2, (j)*b + SEGMENT_SIZE_t5 + 4, v)

#define CBS_R(j) (static_cast<uint32_t>(B[(j)*b + SEGMENT_SIZE_t5 + 4] > B[(j)*b + SEGMENT_SIZE_t5 + 5]))
#define CBS_W(j,v) WriteBlock(B, (j)*b + SEGMENT_SIZE_t5 + 4, (j)*b + SEGMENT_SIZE_t5 + 6, v)

#define TB2_R(j) ReadBlock(B, (j)*b + (bdiv2), (j)*b + (bdiv2) +  SEGMENT_SIZE)
#define TB2_W(j,v) WriteBlock(B, (j)*b + (bdiv2), (j)*b + (bdiv2) +  SEGMENT_SIZE, v)

#define EB2_R(j) (static_cast<uint32_t>(B[(j)*b + (bdiv2) + SEGMENT_SIZE_t3 + 4] > B[(j)*b + SEGMENT_SIZE_t3 + (bdiv2) + 5]))
#define EB2_W(j,v) WriteBlock(B, (j)*b + (bdiv2) + SEGMENT_SIZE_t3 + 4, (j)*b + (bdiv2) + SEGMENT_SIZE_t3 + 6, v)


#define GET_ENDPOINT_A(j) A[(j)*b + (b - 1)]
#define GET_ENDPOINT_B(j) B[(j)*b + (b - 1)]

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


void SetUp(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t b) {
    uint32_t nbA = A.size()/b;
    uint32_t nbB = B.size()/b;

    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        bool in_A = i < nbA;
        uint32_t low = 0;
        uint32_t high = in_A ? nbB : nbA;
        uint32_t k = i - nbA;
        uint32_t mid;

        while (low < high) {
            mid = low + (high-low)/2;
            if ((in_A && GET_ENDPOINT_A(i) > GET_ENDPOINT_B(mid)) 
                 || (!in_A && GET_ENDPOINT_B(k) > GET_ENDPOINT_A(mid))) low = mid+1;
            else high = mid;
        }

        if (in_A) RA_W(i, low);
        else RB_W(k, low);
    });

    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        uint32_t k = i - nbA;
        if (i < nbA) TA_W(i, RA_R(i) + i);
        else TB_W(k, RB_R(k) + k);
    });

    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        if (i < nbA) {
            uint32_t r = RA_R(i);
            uint32_t v = (r >= nbB) ? nbA + nbB : TB_R(r);
            invA_W(i, v);
        } else {
            uint32_t k = i - nbA;
            uint32_t r = RB_R(k);
            uint32_t v = (r >= nbA) ? nbA + nbB : TA_R(r);
            invB_W(k, v);
        }
    });
}

inline bool Done(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t b, bool* flag) {
    std::atomic<bool> done(true);
    uint32_t nbA = A.size()/b;
    uint32_t nbB = B.size()/b;
    parlay::parallel_for(0, nbA + nbB, [&](uint32_t i) {
        if (i < nbA && DA_R(i) == 0 || i >= nbA && DB_R(i - nbA) == 0) done.store(false, std::memory_order_relaxed);
    });
    *flag = done.load(std::memory_order_relaxed);
    
    //*flag = true;
    // parlay::parallel_for(0, seq.size()/b, [&] (uint32_t i) {
    //     if (D_R(i) == 0) *flag = false;
    // });

    return *flag;
}

void EndSort(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t b, bool* flag) {
    uint32_t nbA = A.size()/b;
    uint32_t nbB = B.size()/b;

    std::atomic<bool> done(true);
    std::random_device rd;
    std::size_t seed = (static_cast<std::size_t>(rd()) << 32) ^ rd();
    parlay::random_generator gen(seed);
    std::uniform_int_distribution<uint64_t> dis(0, std::numeric_limits<uint64_t>::max());
    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        auto r = gen[i];
        uint64_t x = dis(r);
        if (i < nbA) {
            CA_W(i, x);
            EA_W(i, 0);
            if (TA_R(i) == i) DA_W(i, 1);
            else {
                DA_W(i, 0);
                done.store(false, std::memory_order_relaxed);
            } 
        } else {
            uint32_t k = i - nbA;
            EB_W(k, 0);
            CB_W(k, x);
            if (TB_R(k) == i) DB_W(k, 1);
            else {
                DB_W(k, 0);
                done.store(false, std::memory_order_relaxed);
            }
        }
    });

    int itcount = 0;
    while (!done.load(std::memory_order_relaxed)) {
        done.store(true, std::memory_order_relaxed);
        if (itcount >= 64) {
            parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
                if (i < nbA) {
                    if (DA_R(i) == 0) {
                        auto r = parlay::random(131542391u 
                                                + itcount * 0x9e3779b97f4a7c15ULL 
                                                + i);
                        uint8_t bit = static_cast<uint8_t>(r.ith_rand(0) & 1ULL);
                        CAS_W(i, bit);
                    }
                } else {
                    uint32_t k = i - nbA;
                    if (DB_R(k) == 0) {
                        auto r = parlay::random(131542391u 
                                                + itcount * 0x9e3779b97f4a7c15ULL 
                                                + i);
                        uint8_t bit = static_cast<uint8_t>(r.ith_rand(0) & 1ULL);
                        CBS_W(k, bit);
                    }
                }
            });
        }

        parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
            if (i < nbA) {
                if (DA_R(i) == 0) { 
                    done.store(false, std::memory_order_relaxed);
                    auto ti = TA_R(i);
                    auto ic = (itcount < 64) ? CA_R(i, itcount) : CAS_R(i);
                    //auto ic = CA_R(i, itcount);
                    if (ic == 1) {
                        uint8_t tc = ti < nbA ? ( itcount < 64 ? CA_R(ti, itcount) : CAS_R(ti)) : ( itcount < 64 ? CB_R(ti - nbA, itcount) : CBS_R(ti - nbA));
                        if (tc == 0) {
                            EA_W(i, 1);
                        }
                    }
                    TA2_W(i, ti);
                    EA2_W(i, EA_R(i));
                }
            } else {
                uint32_t k = i - nbA;
                if (DB_R(k) == 0) {
                    done.store(false, std::memory_order_relaxed);
                    auto tk = TB_R(k);
                    // auto ic = (itcount < 64) ? CB_R(k, itcount) : CBS_R(k);
                    auto ic = (itcount < 64) ? CB_R(k, itcount) : CBS_R(k);
                    if (ic == 1) {
                        uint8_t tc = tk < nbA ? ( itcount < 64 ? CA_R(tk, itcount) : CAS_R(tk)) : ( itcount < 64 ? CB_R(tk - nbA, itcount) : CBS_R(tk - nbA));
                        //uint8_t tc = tk < nbA ? CA_R(tk, itcount) : CB_R(tk - nbA, itcount);
                        if (tc == 0) {
                            EB_W(k, 1);
                        }
                    }
                    TB2_W(k, tk);
                    EB2_W(k, EB_R(k));
                }
            }
        });

        parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
            if (i < nbA) {
                if (EA2_R(i) == 1) {
                    DA_W(i, 1);
                    RA_W(i, i);
                    auto t = TA2_R(i);
                    uint32_t k1 = (t < nbA) ? t : t - nbA;
                    auto& D = (t < nbA) ? A : B;
                    SwapBlockCpy(A, D, i*b, i*b + bdiv2, k1*b, k1*b + bdiv2);
                }
            } else {
                uint32_t k = i - nbA;
                if (EB2_R(k) == 1) {
                    DB_W(k, 1);
                    RB_W(k, i);
                    auto t = TB2_R(k);
                    uint32_t k1 = (t < nbA) ? t : t - nbA;
                    auto& D = (t < nbA) ? A : B;
                    SwapBlockCpy(B, D, k*b, k*b + bdiv2, k1*b, k1*b + bdiv2);
                }
            }
        });

        parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
            if (i < nbA) {
                if (EA_R(i) == 1) {
                    EA_W(i, 0);
                    auto t = RA_R(i);
                    if (t < nbA) {
                        if (TA_R(t) == t) DA_W(t, 1);
                    } else {
                        if (TB_R(t - nbA) == t) DB_W(t - nbA, 1);
                    }
                    uint32_t k1 = (t < nbA) ? t : t - nbA;
                    auto& D = (t < nbA) ? A : B;
                    SwapBlockCpy(A, D, i*b + bdiv2, i*b + b, k1*b + bdiv2, k1*b + b);
                }
            } else {
                uint32_t k = i - nbA;
                if (EB_R(k) == 1) {
                    EB_W(k, 0);
                    auto t = RB_R(k);
                    if (t < nbA) {
                        if (TA_R(t) == t) DA_W(t, 1);
                    } else {
                        if (TB_R(t - nbA) == t) DB_W(t - nbA, 1);
                    }
                    uint32_t k1 = (t < nbA) ? t : t - nbA;
                    auto& D = (t < nbA) ? A : B;
                    SwapBlockCpy(B, D, k*b + bdiv2, k*b + b, k1*b + bdiv2, k1*b + b);
                }
            }
     
        });
        itcount++;
    }
}

// Assumes nb >= 2
inline void Separate(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t start, uint32_t end, uint32_t b, bool base_case) { 
    uint32_t nBlocks = (end - start) / b;
    uint32_t nbA = A.size() / b;
    uint32_t i = (start / b) + (nBlocks / 2) - 1;
    uint32_t e = start/b + nBlocks - 1;
    uint32_t Asize = A.size();
    uint32_t i_index = i*b;
    uint32_t i_inv = i_index < Asize ? invA_R(i) : invB_R(i - nbA);
    if (i_inv > e && !base_case) return;
    uint32_t j = std::min(i_inv, e);
    uint32_t j_index = j*b;

    uint32_t inv1 = i_inv;
    uint32_t inv2 = j_index < Asize ? invA_R(j) : invB_R(j - nbA);
    auto& D = i_index < Asize ? A : B;
    auto seq_index = i_index < Asize ? i_index : i_index - Asize;
    auto D1 = parlay::make_slice(D.begin() + seq_index, D.begin() + seq_index + b);
    auto& D0 = j_index < Asize ? A : B;
    seq_index = j_index < Asize ? j_index : j_index - Asize;
    auto D2 = parlay::make_slice(D0.begin() + seq_index, D0.begin() + seq_index + b);

    if (b <= 128) {
        BubbleSort(D1, D2);  
    } else {
        PairwiseSort(D1);
        PairwiseSort(D2);
        merge(D1, D2);
    }

    if (!base_case) {
        if (i_index < Asize) invA_W(i, inv1);
        else invB_W(i - nbA, inv1);

        if (j_index < Asize) invA_W(j, inv2);
        else invB_W(j - nbA, inv2);
    }
}

void SeqSort(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t start, uint32_t end, uint32_t b) {
    uint32_t n = end - start;
    uint32_t nBlocks = n / b;   
    uint32_t halfBlocks  = nBlocks / 2;
    uint32_t mid = start + (halfBlocks * b); 

    if (n == 2*b) {
        Separate(A, B, start, end, b, true);
        return;
    } else if (n == b) {
        uint32_t Asize = A.size();
        auto& D = end < Asize ? A : B;
        auto x = end < Asize ? 0 : Asize;
        PairwiseSort(parlay::make_slice(D.begin() + start - x, D.begin() + end - x));
        return;
    }
    
    Separate(A, B, start, end, b, false);

    parlay::par_do( 
        [&] {
            SeqSort(A, B, start, mid, b);
        },
        [&] {
            SeqSort(A, B, mid, end, b);
        }
    );
}

bool CheckInversionPointers(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t b) {
    uint32_t nbA = A.size()/b;
    uint32_t nbB = B.size()/b;

    for (int i = 0; i < nbA; i++) {
        auto e_i = GET_ENDPOINT_A(i);
        auto inv_i = RA_R(i);

        for (int j = 0; j < nbB; j++) {
            auto e_j = GET_ENDPOINT_B(j);
            if (j < inv_i) {
                if (e_i < e_j) return false;
                
            } else {
                if (e_i > e_j) return false;
            }
        }  
    }

    for (int i = 0; i < nbB; i++) {
        auto e_i = GET_ENDPOINT_B(i);
        auto inv_i = RB_R(i);

        for (int j = 0; j < nbA; j++) {
            auto e_j = GET_ENDPOINT_A(j);
            if (j < inv_i) {
                if (e_i < e_j) return false;
                
            } else {
                if (e_i > e_j) return false;
            }
        }  
    }

    return true;
}

bool CheckEndSorted(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t b) {
    uint32_t nbA = A.size()/b;
    uint32_t nbB = B.size()/b;

    for (int i = 0; i < nbA + nbB - 1; i++) {
        if (i < nbA - 1) {
            if (GET_ENDPOINT_A(i) > GET_ENDPOINT_A(i+1)) return false;
        } else if (i == nbA - 1) {
            if (GET_ENDPOINT_A(nbA - 1) > GET_ENDPOINT_B(0)) return false;
        } else {
            if (GET_ENDPOINT_B(i - nbA) > GET_ENDPOINT_B(i+1 - nbA)) return false;
        }
    }

    return true;
}

void Merge(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t b) {
    assert(A.size() >= b);
    assert(B.size() >= b);
    assert(A.size() % b == 0);
    assert(B.size() % b == 0);
    assert(A.size() >= b);
    assert(b % 2 == 0);
    assert(b >= 5*SEGMENT_SIZE);

    bdiv2 = b/2;

    bool* flag = (bool*) std::malloc(sizeof(bool));
    *flag = false;
    // std::cout << "Setting Up...\n\n";
    //  auto start = std::chrono::high_resolution_clock().now();
    SetUp(A, B, b);
    //  auto end = std::chrono::high_resolution_clock().now();
    // auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    //assert(CheckInversionPointers(A, B, b));
    // std::cout << "Time for SetUp: " << time.count() << " \n";
    // // std::cout << "End Sorting...\n\n";
    //start = std::chrono::high_resolution_clock().now();
    EndSort(A, B, b, flag);
    // end = std::chrono::high_resolution_clock().now();
    // time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    //assert(CheckEndSorted(A, B, b));
    //std::cout << "Time for EndSort: " << time.count() << " \n";
    //std::cout << "Seq Sorting...\n\n";
    //start = std::chrono::high_resolution_clock().now();
    SeqSort(A, B, 0, (uint32_t)A.size() + B.size(), b);
    // end = std::chrono::high_resolution_clock().now();
    // time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    // std::cout << "Time for SeqSort: " << time.count() << " \n";
    std::free(flag);
}
