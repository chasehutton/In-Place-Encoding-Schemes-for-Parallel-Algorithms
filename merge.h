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

uint32_t A_mod_b, B_mod_b;
uint32_t left_phantom_inversion_pointer;
uint32_t left_phantom_rank;
uint32_t left_phantom_end_sorted_position;
uint32_t first_by_endpoint;
uint32_t first_block_inversion_pointer;
int sub;

void SetUp(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
           parlay::sequence<uint32_t>& left_random, parlay::sequence<uint32_t>& right_random,
           uint32_t tainted_block_index) {
    uint32_t nbA = A.size()/b + (A_mod_b == 0 ? 0 : 1);
    uint32_t nbB = B.size()/b  + (B_mod_b == 0 ? 0 : 1);
    auto x = B_mod_b == 0 ? 0 : 1;

    parlay::parallel_for(0, nbA + nbB - x, [&] (uint32_t i) {
        bool in_A = i < nbA;
        uint32_t low = 0;
        uint32_t high = in_A ? nbB - (B_mod_b == 0 ? 0 : 1) : nbA;
        uint32_t k = i - nbA;
        uint32_t mid;
        uint32_t left_phantom_endpoint = (A_mod_b != 0) ? A[A_mod_b - 1] : get_endpoint(A, 0, sub);
        uint32_t e_a = in_A ? (i == 0 ? left_phantom_endpoint : get_endpoint(A, i, sub)) : 0;
        uint32_t e_b = !in_A ? get_endpoint(B, k, 0) : 0;

        // potential issue is that high might be right_phantom
        while (low < high) {
            mid = low + (high-low)/2;
            e_a = (!in_A) ? (mid == 0 ? left_phantom_endpoint : get_endpoint(A, mid, sub)) : e_a;
            if ((in_A && e_a > get_endpoint(B, mid, 0)) 
                 || (!in_A && e_b > e_a)) low = mid+1;
            else high = mid;
        }

        if (in_A) {
            if (i == 0) left_phantom_rank = low;
            else write_rank(A, i, low, sub);
        }
        else write_rank(B, k, low, 0);
    });

    left_phantom_end_sorted_position = left_phantom_rank;

    first_by_endpoint = 0;

    parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
        uint32_t k = i - nbA; 
        if (i < nbA) {
            auto y = read_rank(A, i, sub) + i;
            if (x == 0) first_by_endpoint = i;
            setup_end_sorted_position_buffer(A, i, y, sub);
        } 
        else {
            auto y =  read_rank(B, k, 0) + k;
            if (x == 0) first_by_endpoint = i;
            setup_end_sorted_position_buffer(B, k, y, 0);
        } 
    });

    parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
        if (i < nbA) {
            auto y = B_mod_b == 0 ? nbB : nbB - 1;
            uint32_t r = read_rank(A, i, sub);
            uint32_t v = (r == x) ? nbA + y : B[r*b + bdiv2 + PADDING];
            write_inversion_pointer(A, i, v, sub);
        } else {
            uint32_t k = i - nbA;
            uint32_t r = read_rank(B, k, 0);
            uint32_t v = (r >= nbA) ? nbA + nbB : (r == 0 ? left_phantom_end_sorted_position : A[r*b + bdiv2 + PADDING + sub]);
            write_inversion_pointer(B, k, v, 0);
        }
    });

    left_phantom_inversion_pointer = B[left_phantom_end_sorted_position*b + bdiv2 + PADDING];
}

void EndSort(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
             parlay::sequence<uint32_t>& left_phantom, parlay::sequence<uint32_t>& right_phantom,
             uint32_t* tainted_block_index, uint32_t num_iterations) {
    assert(num_iterations <= 32);
    uint32_t nbA = A.size()/b  + (A.size() % b == 0 ? 0 : 1);
    uint32_t nbB = B.size()/b  + (B.size() % b == 0 ? 0 : 1);
    auto x = B_mod_b == 0 ? 0 : 1;


    std::random_device rd;
    std::size_t seed = (static_cast<std::size_t>(rd()) << 32) ^ rd();
    parlay::random_generator gen(seed);
    std::uniform_int_distribution<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
    parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
        auto r = gen[i];
        uint32_t random_string = dis(r);
        if (i < nbA) {
            write_coin_flips(A, i, random_string, sub);
            write_swap_flag(A, i, 0, sub);
        } else {
            uint32_t k = i - nbA;
            write_coin_flips(B, k, random_string, 0);
            write_swap_flag(B, k, 0, 0);
        }
    });

    // swap smallest by endpoint with first block
    if (first_by_endpoint != 0) {
        *tainted_block_index = first_by_endpoint;
        auto& S = first_by_endpoint < nbA ? A : B;
        auto st = first_by_endpoint < nbA ? (first_by_endpoint)*b + sub : (first_by_endpoint - nbA)*b;
        auto x = first_by_endpoint < nbA ? sub : 0; 
        auto inv = read_inversion_pointer(S, st, sub);
        first_block_inversion_pointer = left_phantom_inversion_pointer;
        reset_encoding_minus_inv(S, st, x);
        reset_inv_encoding(S, st, x);
        restore_end_sorted_position_buffer(S, st, x);
        swap_block_cpy(A, S, left_phantom, 0, A_mod_b, st, st + b);
        left_phantom_inversion_pointer = inv;
    }

    // swap phantom block with desired position
    if (left_phantom_end_sorted_position != 0) {
        auto& S1 = *tainted_block_index < nbA ? A : B;
        auto st1 = *tainted_block_index < nbA ? *tainted_block_index == 0 ? 0 : (*tainted_block_index)*b + sub : (*tainted_block_index - nbA)*b;
        auto& S2 = left_phantom_end_sorted_position < nbA ? A : B;
        auto x = left_phantom_end_sorted_position < nbA ? sub : 0;
        auto st2 = left_phantom_end_sorted_position < nbA ? (left_phantom_end_sorted_position)*b + sub : (left_phantom_end_sorted_position - nbA)*b;
        swap_block_cpy(S1, S2, st1, st1 + b, st2, st2 + b);
        *tainted_block_index = left_phantom_end_sorted_position;
     }

     std::cout << "Made it 1";
     std::cout.flush();


    for (int it = 0; it < num_iterations; it++) {
        parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
            if (i != *tainted_block_index) {
                if (i < nbA) {
                    auto ti = A[i*b + sub + bdiv2 + PADDING];
                    if (ti != i) {
                        if (read_coin_flip(A, i, it, sub) == 1) {
                            auto tc = ti < nbA ? read_coin_flip(A, ti, it, sub) : read_coin_flip(B, ti - nbA, it, 0);
                            if (tc == 0) {
                                write_swap_flag(A, i, 1, sub);
                                write_rank(A, i, i, sub);
                                auto& D = ti < nbA ? A : B;
                                ti = ti < nbA ? (ti)*b + sub : (ti - nbA)*b;
                                swap_block_cpy_half(A, D, (i)*b + sub, (i)*b + sub + bdiv2, ti, ti + bdiv2);
                            }
                        }
                    }
                } else {
                    uint32_t k = i - nbA;
                    auto tk = B[k*b + bdiv2 + PADDING];
                    if (tk != i) {
                        if (read_coin_flip(B, k, it, 0) == 1) {
                            auto tc = tk < nbA ? read_coin_flip(A, tk, it, sub) : read_coin_flip(B, tk - nbA, it, 0);
                            if (tc == 0) {
                                write_swap_flag(B, k, 1, 0);
                                write_rank(B, k, i, 0);
                                auto& D = tk < nbA ? A : B;
                                tk = tk < nbA ? tk*b + sub : (tk - nbA)*b;
                                swap_block_cpy_half(B, D, k*b, k*b + bdiv2, tk, tk + bdiv2);
                            }
                        }
                    }
                }
            }
        });

        std::cout << "Made it 2";
        std::cout.flush();

        parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
            if (i != *tainted_block_index) {
                if (i < nbA) {
                    if (read_swap_flag(A, i, sub) == 1) {
                        write_swap_flag(A, i, 0, sub);
                        auto t = read_rank(A, i, sub);
                        assert(t != 0);
                        auto& D = t < nbA ? A : B;
                        t = t < nbA ? (t)*b + sub : (t - nbA)*b;
                        swap_block_cpy_half(A, D, (i)*b + sub + bdiv2, (i)*b + sub + b, t + bdiv2, t + b);
                    }

                } else {
                    uint32_t k = i - nbA;
                    if (read_swap_flag(B, k, 0) == 1) {
                        write_swap_flag(B, k, 0, 0);
                        auto t = read_rank(B, k, 0);
                        assert(t != 0);
                        auto& D = t < nbA ? A : B;
                        t -= (t < nbA ? 0 : nbA); 
                        swap_block_cpy_half(B, D, k*b + bdiv2, k*b + b, t*b + bdiv2, t*b + b);
                    }
                }
            }
        });
    } 
    std::cout << "Made it 3";
    std::cout.flush();


    parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
        if (i != *tainted_block_index) {
            if (i < nbA) {
                auto next = A[(i)*b + sub + bdiv2 + PADDING];
                while (next != i) {
                    if (next < i) {
                        mark_self(A, i, 1, sub);
                        break;
                    }
                    next = next < nbA ? A[(next)*b + sub + bdiv2 + PADDING] : B[(next - nbA)*b + bdiv2 + PADDING];
                }
            } else {
                uint32_t k = i - nbA;
                auto next = B[k*b + bdiv2 + PADDING];
                while (next != i) {
                    if (next < i) {
                        mark_self(B, k, 1, 0);
                        break;
                    }
                    next = next < nbA ? A[(next)*b + sub + bdiv2 + PADDING] : B[(next - nbA)*b + bdiv2 + PADDING];
                }
            }
        }
    });

    std::cout << "Made it 4";
    std::cout.flush();


    parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
        if (i != *tainted_block_index) {
            if (i < nbA) {
                auto t = A[(i)*b + sub + bdiv2 + PADDING];
                if (t != i && read_mark(A, i, sub) == 0) {
                    while (t != i) {
                        uint32_t k1 = t < nbA ? (t)*b + sub : (t - nbA)*b;
                        auto& D = t < nbA ? A : B;
                        swap_block_cpy(A, D, (i)*b + sub, (i)*b + sub + b, k1, k1 + b);
                        t = A[(i)*b + sub + bdiv2 + PADDING];
                    }
                }
            } else {
                uint32_t k = i - nbA;
                auto t = B[k*b + bdiv2 + PADDING];
                if (t != i && read_mark(B, k, 0) == 0) {
                    while (t != i) {
                        uint32_t k1 = t < nbA ? (t)*b + sub : (t - nbA)*b;
                        auto& D = t < nbA ? A : B;
                        swap_block_cpy(B, D, k*b, k*b + b, k1, k1 + b);
                        t = B[k*b + bdiv2 + PADDING];
                    }
                }
            }
        }
    });
   
    ///////// comeback to
    parlay::parallel_for(1, nbA + nbB - x, [&] (uint32_t i) {
        if (i != *tainted_block_index) {
            if (i < nbA) {
                restore_end_sorted_position_buffer(A, i, sub);
                reset_encoding_minus_inv(A, i, sub);
            } else {
                uint32_t k = i - nbA;
                restore_end_sorted_position_buffer(B, k, 0);
                reset_encoding_minus_inv(B, k, 0);
            }
        }
    });
}


// inline void Separate(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, parlay::sequence<uint32_t>& right_phantom, uint32_t* tainted_block_index,
//                      uint32_t start, uint32_t end, bool base_case) { 
//     uint32_t nBlocks = (end - start) / b;
//     uint32_t nbA = A.size() / b;
//     uint32_t nbB = B.size() / b;  
//     uint32_t i = (start / b) + (nBlocks / 2) - 1;
//     uint32_t e = start/b + nBlocks - 1;
//     uint32_t Asize = A.size();
//     uint32_t i_index = i*b;
//     uint32_t i_inv;
//     uint32_t j_inv;
//     uint32_t j;
//     uint32_t j_index;
//     uint32_t seq_index;
    
//     if (i == *tainted_block_index) {
//         i_inv = left_phantom_inversion_pointer;
//         if (i_inv > e && !base_case) return;
//         j = std::min(i_inv, e);
//         j_index = j*b;
//         // normal j
//         j_inv = j_index < Asize ? read_inversion_pointer(A, j, A_mod_b) : read_inversion_pointer(B, j - nbA, 0);
//     } else { // i is not first block and is regular
//         i_inv = i_index < Asize ? read_inversion_pointer(A, i, A_mod_b) : read_inversion_pointer(B, i - nbA, 0);
//         j = std::min(i_inv, e);
//         j_index = j*b;
//         // i is not first block and is regular and targets prefix tainted block
//         if (j == *tainted_block_index) {
//             if (i_inv > e && !base_case) return;
//             j_inv = left_phantom_inversion_pointer;   
//         } else if (j == nbA + nbB) { // i targets suffix tainted block
//     } else { // regular case
//         if (i_inv > e && !base_case) return;
//         j_inv = j_index < Asize ? read_inversion_pointer(A, j, A_mod_b) : read_inversion_pointer(B, j - nbA, 0);
//         }
//     }

   
//     // pairwise_sort(D1);
//     // pairwise_sort(D2);

//     if (inv2 != 0) {    
//         parlay::par_do( 
//             [&] {
//                 reset_inv_encoding(D1);
//             },
//             [&] {
//                 reset_inv_encoding(D2);
//             }
//         );
//     } else {
//         reset_inv_encoding(D1);
//     }

//     if (D1[b-1] > D2[0]) merge(D1, D2);
  

//     if (!base_case) {
//         // if (i_index < Asize) write_inversion_pointer(A, i, inv1);
//         // else write_inversion_pointer(B, i - nbA, inv1);

//         if (i_index < Asize) mark_self(A, i, 1);
//         else mark_self(B, i - nbA, 1);

//         if (inv2 != 0) {
//             if (j_index < Asize) write_inversion_pointer(A, j, inv2);
//             else write_inversion_pointer(B, j - nbA, inv2);
//         }
//     }
// }



// inline void Separate(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
//                      parlay::sequence<uint32_t>& left_phantom, parlay::sequence<uint32_t>& right_phantom,
//                      uint32_t* tainted_block_index, uint32_t start, uint32_t end, bool base_case) { 
//     uint32_t nBlocks = (end - start) / b;
//     uint32_t nbA = A.size() / b;
//     uint32_t nbB = B.size() / b;
//     uint32_t i;
//     if (start >= A.size()) {
//         i = (start / b) + (nBlocks / 2) - 1;
//     } else {

//     }
//     auto h = start >= A.size() ?    
//     uint32_t i = (start / b) + (nBlocks / 2) - 1;
//     uint32_t e = start/b + nBlocks - 1;
//     if (end == A.size() + B.size() && B.size() % b != 0) e += 1;
//     uint32_t Asize = A.size();
//     uint32_t i_index = i*b;
//     uint32_t i_inv;
//     uint32_t j_inv;
//     uint32_t j;
//     uint32_t j_index;
//     uint32_t seq_index;
//     // i is first block
//     if (i == 0) {
//         i_inv = first_block_inversion_pointer;
//         if (i_inv > e && !base_case) return;
//         j = std::min(i_inv, e);
//         j_index = j*b;
//         auto& D0 = j_index < Asize ? A : B;
//         seq_index = j_index < Asize ? j_index : j_index - Asize;
//         auto D2 = parlay::make_slice(D0.begin() + seq_index, D0.begin() + seq_index + b);

//         // i is first block and tainted
//         if (i == *tainted_block_index) {
//             // normal j
//             j_inv = j_index < Asize ? read_inversion_pointer(A, j, A_mod_b) : read_inversion_pointer(B, j - nbA, 0);
//             // reset j's encoding
//             reset_inv_encoding(D0, seq_index, j_index < Asize ? A_mod_b : 0);
//             merge(parlay::make_slice(A.begin(), A.begin() + A_mod_b), D0);
//             if (!base_case) {
//                 if (j_index < Asize) write_inversion_pointer(A, j, j_inv, A_mod_b);
//                 else write_inversion_pointer(B, j - nbA, j_inv, 0);
//             }
//         } else { // i is first block but not tainted
//             // i targets tainted block
//             if (j == *tainted_block_index) {
//                 j_inv = left_phantom_inversion_pointer;
//                 merge()

//             } else { // i targets regular block
//                 // normal j
//                 j_inv = j_index < Asize ? read_inversion_pointer(A, j, A_mod_b) : read_inversion_pointer(B, j - nbA, 0);
//             }
//         }

//     } else { // i is not first block
//         // i is not first block and is tainted
//         if (i == *tainted_block_index) {
//             i_inv = left_phantom_inversion_pointer;
//             if (i_inv > e && !base_case) return;
//             j = std::min(i_inv, e);
//             j_index = j*b;
//             // normal j
//             j_inv = j_index < Asize ? read_inversion_pointer(A, j, A_mod_b) : read_inversion_pointer(B, j - nbA, 0);
            

//         } else { // i is not first block and is regular
//             i_inv = i_index < Asize ? read_inversion_pointer(A, i, A_mod_b) : read_inversion_pointer(B, i - nbA, 0);
//             j = std::min(i_inv, e);
//             j_index = j*b;
//             // i is not first block and is regular and targets prefix tainted block
//             if (j == *tainted_block_index) {
//                 if (i_inv > e && !base_case) return;
//                 j_inv = left_phantom_inversion_pointer;
                
//             } else if (j == nbA + nbB) { // i targets suffix tainted block

//             } else { // regular case
//                 if (i_inv > e && !base_case) return;
//                 j_inv = j_index < Asize ? read_inversion_pointer(A, j, A_mod_b) : read_inversion_pointer(B, j - nbA, 0);

//             }

//         }
//     }
    

//     auto& D = i_index < Asize ? A : B;
//     auto seq_index = i_index < Asize ? i_index : i_index - Asize;
//     auto D1 = parlay::make_slice(D.begin() + seq_index, D.begin() + seq_index + b);
//     auto& D0 = j_index < Asize ? A : B;
//     seq_index = j_index < Asize ? j_index : j_index - Asize;
//     auto D2 = parlay::make_slice(D0.begin() + seq_index, D0.begin() + seq_index + b);

   
//     pairwise_sort(D1);
//     pairwise_sort(D2);
//     merge(D1, D2);
  

//     if (!base_case) {
//         if (i_index < Asize) write_inversion_pointer(A, i, inv1);
//         else write_inversion_pointer(B, i - nbA, inv1);

//         if (j_index < Asize) write_inversion_pointer(A, j, inv2);
//         else write_inversion_pointer(B, j - nbA, inv2);
//     }
// }

// void SeqSort(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, 
//              parlay::sequence<uint32_t>& right_phantom, uint32_t* tainted_block_index, uint32_t start, uint32_t end) {

//     uint32_t n = end - start;
//     uint32_t nBlocks = n / b;   
//     uint32_t halfBlocks  = nBlocks / 2;
//     uint32_t mid = start + (halfBlocks * b); 

//     if (n == 2*b) {
//         Separate(A, B, right_phantom, tainted_block_index, start, end, true);
//         return;
//     } else if (n == b) {
//         uint32_t Asize = A.size();
//         auto& D = end < Asize ? A : B;
//         auto x = end < Asize ? 0 : Asize;
//         pairwise_sort(parlay::make_slice(D.begin() + start - x, D.begin() + end - x));
//         return;
//     }
    
//     Separate(A, B, right_phantom, tainted_block_index, start, end, false);

//     parlay::par_do( 
//         [&] {
//             SeqSort(A, B, right_phantom, tainted_block_index, start, mid);
//         },
//         [&] {
//             SeqSort(A, B, right_phantom, tainted_block_index, mid, end);
//         }
//     );
// }


bool CheckInversionPointers(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B) {
    uint32_t nbA = A.size()/b + (A_mod_b == 0 ? 0 : 1);
    uint32_t nbB = B.size()/b  + (B_mod_b % b == 0 ? 0 : 1);

    auto x = A_mod_b != 0 ? A[A_mod_b - 1] : get_endpoint(A, 0, 0);
    auto inv_i = left_phantom_rank;

    for (int j = 0; j < nbB; j++) {
        auto e_j = get_endpoint(B, j, 0);
        if (j < inv_i) {
            if (x < e_j) return false;
            
        } else {
            if (x > e_j) return false;
        }
    }  

    for (int i = 1; i < nbA; i++) {
        auto e_i = get_endpoint(A, i, sub);
        auto inv_i = read_rank(A, i, sub);

        for (int j = 0; j < nbB; j++) {
            auto e_j = get_endpoint(B, j, 0);
            if (j < inv_i) {
                if (e_i < e_j) return false;
                
            } else {
                if (e_i > e_j) return false;
            }
        }  
    }

    for (int i = 1; i < nbB; i++) {
        auto e_i = get_endpoint(B, i, 0);
        auto inv_i = read_rank(B, i, 0);

        for (int j = 0; j < nbA; j++) {
            auto e_j = j == 0 && A_mod_b != 0 ? A[A_mod_b - 1] : get_endpoint(A, j, sub);
            if (j < inv_i) {
                if (e_i < e_j) return false;
                
            } else {
                if (e_i > e_j) return false;
            }
        }  
    }

    return true;
}

bool CheckEndSorted(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t tainted_block_index) {
    uint32_t nbA = A.size()/b  + (A.size() % b == 0 ? 0 : 1);
    uint32_t nbB = B.size()/b  + (B.size() % b == 0 ? 0 : 1);

    for (int i = 1; i < nbA + nbB - 2; i++) {
        if (i < nbA - 1) {
            if (get_endpoint(A, i, sub) > get_endpoint(A, i+1, sub)) return false;
        } else if (i == nbA - 1) {
            if (get_endpoint(A, nbA - 1, sub) > get_endpoint(B, 0,0)) return false;
        } else {
            if (get_endpoint(B, i - nbA, 0) > get_endpoint(B, i+1 - nbA, 0)) return false;
        }
    }

    if (B_mod_b == 0) {
        if (B[(nbB - 2)*b + b - 1] > B[(nbB - 1)*b + b - 1]) return false;
    }

    if (A_mod_b == 0) {
        if (A[b-1] > A[2*b - 1]) return false;
    } else {
        if (A[A_mod_b - 1] > A[2*b - 1 + sub]) return false;
    }

    return true;
}

void Merge(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B) {
    uint32_t Asize = A.size();
    uint32_t Bsize = B.size();

    if (Asize < b && Bsize < b) {
        merge(parlay::make_slice(A), parlay::make_slice(B));
        return;
    }

    if (Asize > Bsize) return Merge(B, A);

    assert(Asize > 2 && Bsize >= 2);

    uint32_t num_iterations = parlay::log2_up(Asize + Bsize / b) ;    

    B_mod_b = Bsize % b;
    A_mod_b = Asize % b;

    uint32_t nbA = A.size()/b  + (A.size() % b == 0 ? 0 : 1);
    uint32_t nbB = B.size()/b  + (B.size() % b == 0 ? 0 : 1);

    sub = (A_mod_b != 0) ? A_mod_b - b : 0;

    parlay::sequence<uint32_t> left_phantom;
    parlay::sequence<uint32_t> right_phantom;
    left_phantom.resize(B_mod_b);
    right_phantom.resize(A_mod_b);

    uint32_t tainted_block_index = 0;

    SetUp(A, B, left_phantom, right_phantom, tainted_block_index);
    // assert(CheckInversionPointers(A,B));
    EndSort(A, B, left_phantom, right_phantom, &tainted_block_index, num_iterations);

    // auto offset = left_phantom_end_sorted_position < nbA ? sub : 0;
    // auto& T = left_phantom_end_sorted_position < nbA ? A : B;
    
    // (left_phantom_end_sorted_position)*b + sub : left_phantom_end_sorted_position*b;
    // auto inv = read_inversion_pointer(T, start, offset);

    assert(CheckEndSorted(A, B, tainted_block_index));


    // if (tainted_block_index == 0) {
    //     auto D = parlay::make_slice(T.begin() + start, T.begin() + start + b);
    //     reset_inv_encoding(D);
    //     merge(parlay::make_slice(A.begin(), A.begin() + A_mod_b), D);
    //     write_inversion_pointer(T, start, inv, offset);
    // } else {
    //     if (left_phantom_end_sorted_position == tainted_block_index) {
    //         auto D = parlay::make_slice(T.begin() + start + A_mod_b, T.begin() + start + b);
    //         merge(parlay::make_slice(A.begin(), A.begin() + A_mod_b), D);
    //         merge(parlay::make_slice(left_phantom), D);
    //         tainted_block_index = left_phantom_end_sorted_position;
    //         left_phantom_inversion_pointer = inv;
    //     } else {
    //         auto D = parlay::make_slice(T.begin() + start, T.begin() + start + b);
    //         reset_inv_encoding(D);
    //         merge(parlay::make_slice(left_phantom), D);
    //         merge(parlay::make_slice(A.begin(), A.begin() + A_mod_b), D);
    //         write_inversion_pointer(T, start, inv, offset);
    //     }
    // }

    // SeqSort(A, B, right_phantom, &tainted_block_index, A_mod_b, (uint32_t)A.size() + B.size());    
}