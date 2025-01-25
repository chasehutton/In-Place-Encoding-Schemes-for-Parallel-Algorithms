#pragma once

#include <cstdint>
#include <random>
#include <iostream>
#include <chrono>


#include "utils.h"
#include "buffer.h"

#include "parlay/sequence.h"
#include "parlay/internal/uninitialized_sequence.h"
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
#include "parlay/utilities.h"

#include "parlay/internal/binary_search.h"
#include "parlay/internal/parallel.h"

#define SEGMENT_SIZE 64

// const inline uint32_t log2(uint32_t x) {
// 	uint32_t e = -1;
//  	for(int i = 0; i < 32; i++) {
//  		if (x & 1) e = i;
//  		x >> 1;
//  	}
//  	return e;
// }

// void out_of_place_merge(buffer& R,
// 						parlay::slice<uint32_t*,uint32_t*> A, parlay::slice<uint32_t*,uint32_t*> B) {
// 	uint32_t nA = A.size();
// 	uint32_t nB = B.size();
// 	uint32_t nR = nA+nB;
// 
// 	if (nR < _merge_base) {
// 		seq_merge<assignment_tag>(A, B, R, f);
// 	  }
// 	else if (nA == 0) {
// 	    parallel_for(0, nB, [&](size_t i) {
// 	      R.aux[i] = B[i]);
// 	    });
// 	  }
// 	else if (nB == 0) {
// 	    parallel_for(0, nA, [&](size_t i) {
// 	      R.aux[i] = A[i];
// 	    });
// 	  }
// 	else {
// 	    size_t mA = nA / 2;
// 	    // important for stability that binary search identifies
// 	    // first element in B greater or equal to A[mA]
// 	    size_t mB = binary_search(B, A[mA]);
// 	    if (mB == 0) mA++;  // ensures at least one on each side
// 	    size_t mR = mA + mB;
// 	    auto left = [&]() {
// 	      merge_into<assignment_tag>(A.cut(0, mA), B.cut(0, mB), R.cut(0, mR), f, cons);
// 	    };
// 	    auto right = [&]() {
// 	      merge_into<assignment_tag>(A.cut(mA, nA), B.cut(mB, nB), R.cut(mR, nR), f, cons);
// 	    };
// 	    par_do(left, right, cons);
// 	  }
// }




void Merge(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B) {
	uint32_t n = A.size();
	uint32_t g = 300
    uint32_t r = n/g;

    buffer Buffer = new buffer(
    							r,
								parlay::make_slice(A.begin()-131*r, A.begin()+n-129*r),
								parlay::make_slice(A.begin()-129*r, A.begin()+n-r),
								64);
	Buffer.initialize();

	for (int i = 0; i < 150; i++) {
		auto C = parlay::make_slice(A.begin()+i*r, A.begin()+i*r+r);
		auto D = parlay::make_slice(A.begin()+n-r, A.end());
		auto E = parlay::make_slice(B.begin()+n-i*r-r, B.end()-i*r);
		auto F = parlay::make_slice(B.begin(), B.begin()+r);

		merge_into(D, E, Buffer);
		distribute(Buffer,D,E);
		
		merge_into(C, D, Buffer);
		distribute(Buffer,D,E);
		
		merge_into(D, E, Buffer);
		distribute(Buffer,D,E);
	}
	

    bool* flag = (bool*) std::malloc(sizeof(bool));
    *flag = false;
    std::cout << "Setting Up...\n\n";
    auto start = std::chrono::high_resolution_clock().now();
    SetUp(seq, b);
    auto end = std::chrono::high_resolution_clock().now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "Time for SetUp: " << time.count() << " \n";
    std::cout << "End Sorting...\n\n";
    start = std::chrono::high_resolution_clock().now();
    EndSort(seq, b, flag);
    end = std::chrono::high_resolution_clock().now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "Time for EndSort: " << time.count() << " \n";
    std::cout << "Seq Sorting...\n\n";
    start = std::chrono::high_resolution_clock().now();
    SeqSort(seq, 0, (uint32_t)seq.size(), b);
    end = std::chrono::high_resolution_clock().now();
    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    std::cout << "Time for SeqSort: " << time.count() << " \n";

    std::free(flag);
}

