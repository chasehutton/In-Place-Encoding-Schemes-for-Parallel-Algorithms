#pragma once

#include <cstdint>
#include <random>
#include <iostream>

#include "utils.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"

#define BLOCK_SIZE 64
#define D_R(j) ReadBlock(seq, ((j)*b) + 2*BLOCK_SIZE + 2, ((j)*b) + 2*BLOCK_SIZE + 3)
#define D_W(j,v) WriteBlock(seq, ((j)*b) + 2*BLOCK_SIZE + 2, ((j)*b) + 2*BLOCK_SIZE + 3, v)
#define E_R(j) ReadBlock(seq, ((j)*b) + 2*BLOCK_SIZE + 4, ((j)*b) + 2*BLOCK_SIZE + 5)
#define E_W(j,v) WriteBlock(seq, ((j)*b) + 2*BLOCK_SIZE + 4, ((j)*b) + 2*BLOCK_SIZE + 5, v)
#define C_R(j) ReadBlock(seq, ((j)*b) + 2*BLOCK_SIZE, ((j)*b) + 2*BLOCK_SIZE + 1)
#define C_W(j,v) WriteBlock(seq, ((j)*b) + 2*BLOCK_SIZE, ((j)*b) + 2*BLOCK_SIZE + 1, v)
#define inv_R(j) ReadBlock(seq, ((j)*b) + BLOCK_SIZE, ((j)*b) + 2*BLOCK_SIZE - 1)
#define inv_W(j,v) WriteBlock(seq, ((j)*b) + BLOCK_SIZE, ((j)*b) + 2*BLOCK_SIZE - 1, v)
#define T_R(j) ReadBlock(seq, ((j)*b), ((j)*b) + BLOCK_SIZE - 1)
#define T_W(j,v) WriteBlock(seq, ((j)*b), ((j)*b) + BLOCK_SIZE - 1, v)
#define GET_ENDPOINT(j) seq[(j)*b + (b-1)]

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
    auto nb = seq.size()/(2*b);

    parlay::parallel_for(0, nb, [&] (auto i) {
        uint32_t low = 0;
        uint32_t high = nb;
        uint32_t mid;

        while (low < high) {
            mid = (high+low)/2;
            if (GET_ENDPOINT(i) > GET_ENDPOINT(mid + nb)) low = mid+1;
            else high = mid;
        }
        T_W(i,low);

        low = 0;
        high = nb;

        while (low < high) {
            mid = (high+low)/2;
            if (GET_ENDPOINT(i + nb) > GET_ENDPOINT(mid)) low = mid+1;
            else high = mid;
        }
        T_W(i + nb, low);
    });
 
    parlay::parallel_for(0, nb, [&] (auto i) {
        inv_W(i, T_R(T_R(i)));
        inv_W(i + nb, T_R(T_R(i + nb)));
    });

    parlay::parallel_for(0, nb, [&] (auto i) {
        T_W(i, T_R(i) + i);
        T_W(i + nb, T_R(i + nb) + i);
    });
}

///////////////////////////////////////////////////////////////////////////////
// Determines when EndSort is finished. the size of seq is assumed to be     //
// divisible by and larger than b. b is assumed larger than BLOCK_SIZE.      //
///////////////////////////////////////////////////////////////////////////////

inline bool Done(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {
    *flag = 1;
    parlay::parallel_for(0, seq.size()/b, [&] (auto i) {
        if (D_R(i) == 0) *flag = 0;
    });

    return *flag;
}

void EndSort(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {
    parlay::parallel_for(0, seq.size()/b, [&] (auto i) {
        E_W(i,0); 
        if (T_R(i) == i) D_W(i,1);
        else D_W(i,0);
    });

    int itcount = 0;
    while(!Done(seq, b, flag)) {
        // if (itcount > 1000) {
        //     for (int i = 0; i < seq.size()/b; i++) {
        //         std::cout << "block: " << i << " data: D_R(i) " << D_R(i) << " T_R(i) " << T_R(i) <<  " C_R(i) " << C_R(i) << " E_R(i) " << E_R(i) << " \n";  
        //     }

        //     return;
        // }
        itcount++;
        parlay::parallel_for(0, seq.size()/b, [&] (auto i) {
            if (D_R(i) == 0) {
                auto r = parlay::random(131542391u + itcount * 0x9e3779b97f4a7c15ULL + i);
                uint8_t bit = static_cast<uint8_t>(r.ith_rand(0) & 1ULL);
                C_W(i, bit);
            } 
        });

        parlay::parallel_for(0, seq.size()/b, [&] (auto i) {
            if (D_R(i) == 0 && C_R(i) == 1 && C_R(T_R(i)) == 0) {
                E_W(i,1);
            }
        });

        parlay::parallel_for(0, seq.size()/b, [&] (auto i) {
            if (D_R(i) == 0 && E_R(i) == 1) {
                D_W(i,1);
                auto t = T_R(i);
                SwapBlock(seq, i*b, i*b+b-1, t*b, t*b+b-1);
                if (T_R(i) == i) D_W(i,1);
            }
        });
     }

}

void Separate(parlay::sequence<uint32_t>& seq, uint32_t b) {

}

void SeqSort(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {

}

void Merge(parlay::sequence<uint32_t>& seq, uint32_t b) {
    assert(seq.size() % b == 0);
    assert(seq.size() > b);
    assert((seq.size()/2) % b == 0);
    assert(b % 2 == 0);
    assert(b >= 4*BLOCK_SIZE);


    bool* flag = (bool*) std::malloc(sizeof(bool));
    *flag = 0;

    std::cout << "Setting Up...\n\n";
    SetUp(seq, b);
    std::cout << "End Sorting...\n\n";
    EndSort(seq, b, flag);
    //SeqSort(seq, b, flag);

    std::free(flag);
}