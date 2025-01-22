#pragma once

#include <cstdint>

#include "utils.h"
#include "parlay/sequence.h"
#include "parlay/parallel.h"

#define BLOCK_SIZE 64
#define D_R(i) ReadBlock(seq, (i) + 2*BLOCK_SIZE + 2, (i) + 2*BLOCK_SIZE + 3)
#define D_W(i,v) WriteBlock(seq, (i) + 2*BLOCK_SIZE + 2, (i) + 2*BLOCK_SIZE + 3, v)
#define C_R(i) ReadBlock(seq, (i) + 2*BLOCK_SIZE, (i) + 2*BLOCK_SIZE + 1)
#define C_W(i,v) WriteBlock(seq, (i) + 2*BLOCK_SIZE, (i) + 2*BLOCK_SIZE + 1, v)
#define inv_R(i) ReadBlock(seq, (i) + BLOCK_SIZE, (i) + 2*BLOCK_SIZE - 1)
#define inv_W(i,v) WriteBlock(seq, (i) + BLOCK_SIZE, (i) + 2*BLOCK_SIZE - 1, v)
#define T_R(i) ReadBlock(seq, (i), (i) + BLOCK_SIZE - 1)
#define T_W(i,v) WriteBlock(seq, (i), (i) + BLOCK_SIZE - 1, v)

///////////////////////////////////////////////////////////////////////////////////////////
// Each worker thread gets a workspace of size c*b where b is the given non-encoding     //
// block size and c is a constant to be determined (currently set to be 4).              //
// The first 3*b uint32_t values of each workspace are reserved for the out-of-place     //
// sequential merging done in the Separate function. The remaining values are used       //
// for any miscellaneous auxilary storage.                                               //
///////////////////////////////////////////////////////////////////////////////////////////

static std::vector<uint32_t*> workspaces;

void Init_workspaces(std::size_t n) {
    auto nw = parlay::num_workers();
    workspaces.resize(nw);
    
    for (int i = 0; i < nw; i++) {
        workspaces[i] = (uint32_t*) std::malloc(n * sizeof(uint32_t));
    }
}

void Free_workspaces() {
    for (auto* ptr : workspaces) {
        std::free(ptr);
    }

    workspaces.clear();
}


void Compute_Ranks(parlay::sequence<uint32_t>& seq, uint32_t b) {
    
    parlay::parallel_for(0, seq.size()/b, [&] (auto i) {

    });
}

///////////////////////////////////////////////////////////////////////////////
// Determines when EndSort is finished. the size of seq is assumed to be     //
// divisible by and larger than b. b is assumed larger than BLOCK_SIZE.      //
///////////////////////////////////////////////////////////////////////////////

bool Done(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {
    *flag = 1;
    parlay::parallel_for(0, seq.size()/b, [&] (auto i) {
        if (D_R(i) == 0) *flag = 0;
    });

    return *flag;
}

void EndSort(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {

}

void Separate(parlay::sequence<uint32_t>& seq, uint32_t b) {

}

void SeqSort(parlay::sequence<uint32_t>& seq, uint32_t b, bool* flag) {

}

void Merge(parlay::sequence<uint32_t>& seq, uint32_t b) {
    assert(seq.size() % 2 == 0);
    assert(seq.size() > b);
    assert((seq.size()/2) % b == 0);
    assert(b >= BLOCK_SIZE);

    // Sequential out-of-place merge will require 3*b and add an additional b for saftey
    Init_workspaces(4*b);
    bool* flag = (bool*) std::malloc(sizeof(bool));
    *flag = 0;

    Compute_Ranks(seq, b);
    SeqSort(seq, b, flag);

    std::free(flag);
    Free_workspaces();
}