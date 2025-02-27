#pragma once

#include <cstdint>
#include <random>
#include <iostream>
#include <chrono>


#include "utils.h"
#include "block_size.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/utilities.h"

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
uint32_t A_mod_b, B_mod_b;
uint32_t left_phantom_inversion_pointer;
uint32_t left_phantom_end_sorted_position;
uint32_t left_phantom_mark;
uint32_t first_by_endpoint;

uint32_t inline get_endpoint(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    return S[j*b + b - 1 + offset];
}

inline void setup_end_sorted_position_buffer(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t v, uint32_t offset) {
   uint32_t seq_index = j*b + bdiv2 + PADDING + offset;
   write_block_64(S, seq_index + 2, S[seq_index]);
   S[seq_index] = v;
}

inline void restore_end_sorted_position_buffer(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    uint32_t seq_index = j*b + bdiv2 + PADDING + offset;
    S[seq_index] = read_block_64(S, seq_index + 2);
}

uint32_t inline read_end_sorted_position(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    uint32_t seq_index = b*j + bdiv2 + PADDING + offset;
    return read_block_64(S, seq_index);
}

void inline write_end_sorted_position(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, uint32_t offset) {
    uint32_t seq_index = b*j + bdiv2 + PADDING + offset;
    write_block_64(S, seq_index, value);
}

uint32_t inline read_inversion_pointer(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    uint32_t seq_index = b*j + offset;
    return read_block_64(S, seq_index);
}

void inline write_inversion_pointer(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, uint32_t offset) {
    uint32_t seq_index = b*j + offset;
    write_block_64(S, seq_index, value);
}

uint32_t inline read_rank(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE + PADDING + offset;
    return read_block_64(S, seq_index);
}

void inline write_rank(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, uint32_t offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE + PADDING + offset;
    write_block_64(S, seq_index, value);
}

uint32_t inline read_coin_flip(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t k, uint32_t offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE + 2 + PADDING_t2 + 2*k + offset;
    return static_cast<uint32_t>(S[seq_index] > S[seq_index + 1]);
}

void inline write_coin_flips(parlay::sequence<uint32_t>& S, uint32_t j, uint64_t r, uint32_t offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE + 2 + PADDING_t2 + offset;
    write_block_64(S, seq_index, r);
}


uint32_t inline read_swap_flag(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE_t2 + PADDING_t2 + offset;
    return static_cast<uint32_t>(S[seq_index] > S[seq_index + 1]);
}

void inline write_swap_flag(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t value, uint32_t offset) {
    uint32_t seq_index = b*j + SEGMENT_SIZE_t2 + PADDING_t2 + offset;
    write_block_2(S, seq_index, value);
}


void inline mark_self(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE_t3 + 2 + PADDING_t3 + offset;
    write_block_2(S, seq_index, 1);
}

uint32_t inline read_mark(parlay::sequence<uint32_t>& S, uint32_t j, uint32_t offset) {
    uint32_t seq_index = b*j + bdiv2 + SEGMENT_SIZE_t3 + 2 + PADDING_t3 + offset;
    return static_cast<uint32_t>(S[seq_index] > S[seq_index + 1]); 
}


void SetUp(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
           parlay::sequence<uint32_t>& left_random, parlay::sequence<uint32_t>& right_random,
           uint32_t tainted_block_index) {
    uint32_t nbA = A.size()/b + 1;
    uint32_t nbB = B.size()/b;

    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        bool in_A = i < nbA;
        uint32_t low = 0;
        uint32_t high = in_A ? nbB - 1 : nbA;
        uint32_t k = i - nbA;
        uint32_t mid;
        uint32_t left_phantom_endpoint = A[A_mod_b - 1];
        uint32_t e_a = in_A ? (i == 0 ? left_phantom_endpoint : get_endpoint(A, i, A_mod_b)) : 0;
        uint32_t e_b = !in_A ? get_endpoint(B, k, 0) : 0;

        // potential issue is that high might be right_phantom
        while (low < high) {
            mid = low + (high-low)/2;
            e_a = (!in_A && mid == 0) ? left_phantom_endpoint : get_endpoint(A, mid, A_mod_b);
            if ((in_A && e_a > get_endpoint(B, mid, 0)) 
                 || (!in_A && e_b > e_a)) low = mid+1;
            else high = mid;
        }

        if (in_A) {
            if (i == 0) left_phantom_end_sorted_position = low;
            else write_rank(A, i, low, A_mod_b);
        }
        else write_rank(B, k, low, 0);
    });

    first_by_endpoint = 0;

    parlay::parallel_for(1, nbA + nbB, [&] (uint32_t i) {
        uint32_t k = i - nbA;
        if (i < nbA) {
            auto x = read_rank(A, i, A_mod_b) + i;
            if (x == 0) first_by_endpoint = i;
            setup_end_sorted_position_buffer(A, i, x, A_mod_b);
        } 
        else {
            auto x =  read_rank(B, k, 0) + k;
            if (x == 0) first_by_endpoint == i;
            setup_end_sorted_position_buffer(B, k, x, 0);
        } 
    });

    parlay::parallel_for(1, nbA + nbB, [&] (uint32_t i) {
        if (i < nbA) {
            uint32_t r = read_rank(A, i, A_mod_b);
            uint32_t v = (r >= nbB) ? nbA + nbB + 1 : B[r*b + bdiv2 + PADDING];
            write_inversion_pointer(A, i, v, A_mod_b);
        } else {
            uint32_t k = i - nbA;
            uint32_t r = read_rank(B, k, 0);
            uint32_t v = (r >= nbA) ? nbA + nbB + 1 : (r == 0 ? left_phantom_end_sorted_position : A[r*b + bdiv2 + PADDING + A_mod_b]);
            write_inversion_pointer(B, k, v, 0);
        }
    });

    left_phantom_inversion_pointer = B[left_phantom_end_sorted_position*b + bdiv2 + PADDING];
}

void EndSort(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
             parlay::sequence<uint32_t>& left_random, parlay::sequence<uint32_t>& right_random,
             uint32_t tainted_block_index, uint32_t num_iterations) {
    assert(num_iterations <= 32);
    uint32_t nbA = A.size()/b + 1;
    uint32_t nbB = B.size()/b;

    std::random_device rd;
    std::size_t seed = (static_cast<std::size_t>(rd()) << 32) ^ rd();
    parlay::random_generator gen(seed);
    std::uniform_int_distribution<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
    parlay::parallel_for(1, nbA + nbB, [&] (uint32_t i) {
        auto r = gen[i];
        uint32_t random_string = dis(r);
        if (i < nbA) {
            write_coin_flips(A, i, random_string, A_mod_b);
            write_swap_flag(A, i, 0, A_mod_b);
        } else {
            uint32_t k = i - nbA;
            write_coin_flips(B, k, random_string, 0);
            write_swap_flag(B, k, 0, 0);
        }
    });

    for (int it = 0; it < num_iterations; it++) {
        parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
            if (i < nbA) {
                auto ti = A[i*b + bdiv2 + PADDING];
                if (ti != i) {
                    if (read_coin_flip(A, i, it) == 1) {
                        auto tc = ti < nbA ? read_coin_flip(A, ti, it) : read_coin_flip(B, ti - nbA, it);
                        if (tc == 0) {
                            write_swap_flag(A, i, 1);
                            write_rank(A, i, i);
                            auto& D = ti < nbA ? A : B;
                            ti -= (ti < nbA ? 0 : nbA);
                            swap_block_cpy_half(A, D, i*b, i*b + bdiv2, ti*b, ti*b + bdiv2);
                        }
                    }
                }
            } else {
                uint32_t k = i - nbA;
                auto tk = B[k*b + bdiv2 + PADDING];
                if (tk != i) {
                    if (read_coin_flip(B, k, it) == 1) {
                        auto tc = tk < nbA ? read_coin_flip(A, tk, it) : read_coin_flip(B, tk - nbA, it);
                        if (tc == 0) {
                            write_swap_flag(B, k, 1);
                            write_rank(B, k, i);
                            auto& D = tk < nbA ? A : B;
                            tk -= (tk < nbA ? 0 : nbA);
                            swap_block_cpy_half(B, D, k*b, k*b + bdiv2, tk*b, tk*b + bdiv2);
                        }
                    }
                }
            }
        });

        parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
            if (i < nbA) {
                if (read_swap_flag(A, i) == 1) {
                    write_swap_flag(A, i, 0);
                    auto t = read_rank(A, i);
                    auto& D = t < nbA ? A : B;
                    t -= (t < nbA ? 0 : nbA); 
                    swap_block_cpy_half(A, D, i*b + bdiv2, i*b + b, t*b + bdiv2, t*b + b);
                }

            } else {
                uint32_t k = i - nbA;
                if (read_swap_flag(B, k) == 1) {
                    write_swap_flag(B, k, 0);
                    auto t = read_rank(B, k);
                    auto& D = t < nbA ? A : B;
                    t -= (t < nbA ? 0 : nbA); 
                    swap_block_cpy_half(B, D, k*b + bdiv2, k*b + b, t*b + bdiv2, t*b + b);
                }
            }
        });
    } 

    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        if (i < nbA) {
            auto next = A[i*b + bdiv2 + PADDING];
            while (next != i) {
                if (next < i) {
                    mark_self(A, i);
                    break;
                }
                next = next < nbA ? A[next*b + bdiv2 + PADDING] : B[(next - nbA)*b + bdiv2 + PADDING];
            }
        } else {
            uint32_t k = i - nbA;
            auto next = B[k*b + bdiv2 + PADDING];
            while (next != i) {
                if (next < i) {
                    mark_self(B, k);
                    break;
                }
                next = next < nbA ? A[next*b + bdiv2 + PADDING] : B[(next - nbA)*b + bdiv2 + PADDING];
            }
        }
    });

    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        if (i < nbA) {
            auto t = A[i*b + bdiv2 + PADDING];
            if (t != i && read_mark(A, i) == 0) {
                while (t != i) {
                    uint32_t k1 = t < nbA ? t : t - nbA;
                    auto& D = t < nbA ? A : B;
                    swap_block_cpy(A, D, i*b, i*b + b, k1*b, k1*b + b);
                    t = A[i*b + bdiv2 + PADDING];
                }
            }
        } else {
            uint32_t k = i - nbA;
            auto t = B[k*b + bdiv2 + PADDING];
            if (t != i && read_mark(B, k) == 0) {
                while (t != i) {
                    uint32_t k1 = t < nbA ? t : t - nbA;
                    auto& D = t < nbA ? A : B;
                    swap_block_cpy(B, D, k*b, k*b + b, k1*b, k1*b + b);
                    t = B[k*b + bdiv2 + PADDING];
                }
            }
        }
    });
   
    parlay::parallel_for(0, nbA + nbB, [&] (uint32_t i) {
        if (i < nbA) {
            restore_end_sorted_position_buffer(A, i);
        } else {
            uint32_t k = i - nbA;
            restore_end_sorted_position_buffer(B, k);
        }
    });
}



inline void Separate(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
                     parlay::sequence<uint32_t>& left_phantom, parlay::sequence<uint32_t>& right_phantom,
                     uint32_t tainted_block_index, uint32_t start, uint32_t end, bool base_case) { 

    uint32_t nBlocks = (end - start) / b;
    uint32_t nbA = A.size() / b;
    uint32_t i = (start / b) + (nBlocks / 2) - 1;
    uint32_t e = start/b + nBlocks - 1;
    uint32_t Asize = A.size();
    uint32_t i_index = i*b;
    uint32_t i_inv = i_index < Asize ? read_inversion_pointer(A, i) : read_inversion_pointer(B, i - nbA);
    if (i_inv > e && !base_case) return;
    uint32_t j = std::min(i_inv, e);
    uint32_t j_index = j*b;

    uint32_t inv1 = i_inv;
    uint32_t inv2 = j_index < Asize ? read_inversion_pointer(A, j) : read_inversion_pointer(B, j - nbA);
    auto& D = i_index < Asize ? A : B;
    auto seq_index = i_index < Asize ? i_index : i_index - Asize;
    auto D1 = parlay::make_slice(D.begin() + seq_index, D.begin() + seq_index + b);
    auto& D0 = j_index < Asize ? A : B;
    seq_index = j_index < Asize ? j_index : j_index - Asize;
    auto D2 = parlay::make_slice(D0.begin() + seq_index, D0.begin() + seq_index + b);

   
    pairwise_sort(D1);
    pairwise_sort(D2);
    merge(D1, D2);
  

    if (!base_case) {
        if (i_index < Asize) write_inversion_pointer(A, i, inv1);
        else write_inversion_pointer(B, i - nbA, inv1);

        if (j_index < Asize) write_inversion_pointer(A, j, inv2);
        else write_inversion_pointer(B, j - nbA, inv2);
    }
}

void SeqSort(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
             parlay::sequence<uint32_t>& left_pantom, parlay::sequence<uint32_t>& right_phantom,
             uint32_t tainted_block_index, uint32_t start, uint32_t end) {

    uint32_t n = end - start;
    uint32_t nBlocks = n / b;   
    uint32_t halfBlocks  = nBlocks / 2;
    uint32_t mid = start + (halfBlocks * b); 

    if (n == 2*b) {
        Separate(A, B, left_pantom, right_phantom, tainted_block_index, start, end, true);
        return;
    } else if (n == b) {
        uint32_t Asize = A.size();
        auto& D = end < Asize ? A : B;
        auto x = end < Asize ? 0 : Asize;
        pairwise_sort(parlay::make_slice(D.begin() + start - x, D.begin() + end - x));
        return;
    }
    
    Separate(A, B, left_pantom, right_phantom, tainted_block_index, start, end, false);

    parlay::par_do( 
        [&] {
            SeqSort(A, B, left_pantom, right_phantom, tainted_block_index, start, mid);
        },
        [&] {
            SeqSort(A, B, left_pantom, right_phantom, tainted_block_index, mid, end);
        }
    );
}


bool CheckInversionPointers(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B) {
    uint32_t nbA = A.size()/b;
    uint32_t nbB = B.size()/b;

    for (int i = 0; i < nbA; i++) {
        auto e_i = get_endpoint(A, i);
        auto inv_i = read_rank(A, i);

        for (int j = 0; j < nbB; j++) {
            auto e_j = get_endpoint(B, j);
            if (j < inv_i) {
                if (e_i < e_j) return false;
                
            } else {
                if (e_i > e_j) return false;
            }
        }  
    }

    for (int i = 0; i < nbB; i++) {
        auto e_i = get_endpoint(B, i);
        auto inv_i = read_rank(B, i);

        for (int j = 0; j < nbA; j++) {
            auto e_j = get_endpoint(A, j);
            if (j < inv_i) {
                if (e_i < e_j) return false;
                
            } else {
                if (e_i > e_j) return false;
            }
        }  
    }

    return true;
}

bool CheckEndSorted(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B) {
    uint32_t nbA = A.size()/b;
    uint32_t nbB = B.size()/b;

    for (int i = 0; i < nbA + nbB - 1; i++) {
        if (i < nbA - 1) {
            if (get_endpoint(A, i) > get_endpoint(A, i+1)) return false;
        } else if (i == nbA - 1) {
            if (get_endpoint(A, nbA - 1) > get_endpoint(B, 0)) return false;
        } else {
            if (get_endpoint(B, i - nbA) > get_endpoint(B, i+1 - nbA)) return false;
        }
    }

    return true;
}

void Merge(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B) {
    uint32_t Asize = A.size();
    uint32_t Bsize = B.size();

    if (Asize > Bsize) return Merge(B, A);

    assert(Asize > 2 && Bsize >= 2);

    uint32_t num_iterations = parlay::log2_up(Asize + Bsize / b) ;    

    B_mod_b = Bsize % b;
    A_mod_b = Asize % b;

    parlay::sequence<uint32_t> left_phantom;
    parlay::sequence<uint32_t> right_phantom;
    left_phantom.resize(B_mod_b);
    right_phantom.resize(A_mod_b);

    uint32_t tainted_block_index = 0;

    SetUp(A, B, left_phantom, right_phantom, tainted_block_index);
    EndSort(A, B, left_phantom, right_phantom, tainted_block_index, num_iterations);
    SeqSort(A, B, left_phantom, right_phantom, tainted_block_index, 0, (uint32_t)A.size() + B.size());    
}