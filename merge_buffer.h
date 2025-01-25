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
#include "parlay/parallel.h"

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


void distribute (	parlay::slice<uint32_t*,uint32_t*> aux,
					parlay::slice<uint32_t*,uint32_t*> A,
					parlay::slice<uint32_t*,uint32_t*> B) {
	uint32_t nA = A.size();
	uint32_t nB = B.size();
	parlay::parallel_for(0, nA, [&](uint32_t i) {
    	A[i] = aux[i];
  	});
  	parlay::parallel_for(0, nB, [&](uint32_t i) {
		B[i] = aux[nA+i];
	});
}

void buffer_merge(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B) {
	std::cout << "Buffer Merging...\n\n";
	auto start = std::chrono::high_resolution_clock().now();
	
	uint32_t n = A.size();
	uint32_t g = 512;
    uint32_t r = n/g;

    buffer Buffer = buffer(
    							r,
								parlay::make_slice(A.begin()-131*r, A.begin()+n-129*r),
								parlay::make_slice(A.begin()-129*r, A.begin()+n-r),
								64);
	Buffer.initialize();

	for (int i = 0; i < g/2; i++) {
		auto C = parlay::make_slice(A.begin()+i*r, A.begin()+i*r+r);
		auto D = parlay::make_slice(A.begin()+n-r, A.end());
		auto E = parlay::make_slice(B.begin(), B.begin()+r);
		auto F = parlay::make_slice(B.begin()+n-i*r-r, B.end()-i*r);

		parlay::internal::merge_into<parlay::move_assign_tag>(D, E, Buffer.aux, std::less<>());
		distribute(Buffer.aux, D, E);
		
		parlay::internal::merge_into<parlay::move_assign_tag>(C, D, Buffer.aux, std::less<>());
		distribute(Buffer.aux, C, D);
		
		parlay::internal::merge_into<parlay::move_assign_tag>(E, F, Buffer.aux, std::less<>());
		distribute(Buffer.aux, E, F);
	}

	Buffer.restore();

	Buffer.aux = parlay::make_slice(A.begin(), A.begin()+2*r);
	Buffer.enc = parlay::make_slice(A.begin()+2*r, A.begin()+130*r),

	Buffer.initialize();

	for (int i = g/2; i < g; i++) {
		auto C = parlay::make_slice(A.begin()+i*r, A.begin()+i*r+r);
		auto D = parlay::make_slice(A.begin()+n-r, A.end());
		auto E = parlay::make_slice(B.begin(), B.begin()+r);
		auto F = parlay::make_slice(B.begin()+n-i*r-r, B.end()-i*r);

		parlay::internal::merge_into<parlay::move_assign_tag>(D, E, Buffer.aux, std::less<>());
		distribute(Buffer.aux, D, E);
		
		parlay::internal::merge_into<parlay::move_assign_tag>(C, D, Buffer.aux, std::less<>());
		distribute(Buffer.aux, C, D);
		
		parlay::internal::merge_into<parlay::move_assign_tag>(E, F, Buffer.aux, std::less<>());
		distribute(Buffer.aux, E, F);
	}

	auto D = parlay::make_slice(A.begin()+n-r, A.end());
	auto E = parlay::make_slice(B.begin(), B.begin()+r);
	parlay::internal::merge_into<parlay::uninitialized_move_tag>(D, E, Buffer.aux, std::less<>());
	distribute(Buffer.aux, D, E);

	Buffer.restore();
	
    auto end = std::chrono::high_resolution_clock().now();
    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);
    
    std::cout << "Time for Buffer Merge: " << time.count() << " \n";
}

