#pragma once

#include <cstdint>
#include <cstring> 


#include "parlay/sequence.h"
#include "parlay/parallel.h"
#include "utils.h"

struct buffer {
	uint32_t aux_size;
	uint32_t segment_size;
	parlay::slice<uint32_t*,uint32_t*> aux;
	parlay::slice<uint32_t*,uint32_t*> enc;

	 // Constructor
    buffer(uint32_t aux_size, parlay::slice<uint32_t*,uint32_t*> aux, parlay::slice<uint32_t*,uint32_t*> enc, uint32_t segment_size) 
        : aux_size(aux_size), 
          segment_size(segment_size),
          aux(aux), 
          enc(enc) {
          // check size parameters match up
          // assert(segment_size == 64 && "incorrect segment size");
          // assert(aux.size() == aux_size && "aux size mismatch");
          // assert(enc.size() == aux_size()*segment_size && "enc size mismatch");
    }

    uint32_t read(uint32_t i) {
    	// assert(0 <=  i && i < aux_size()) && "aux out of bounds error";
  		return ReadBlock(enc, i*segment_size, (i+1)*segment_size);
    }
    
    void write(uint32_t i, uint32_t v) {
    	// assert(0 <=  i && i < aux_size()) && "aux out of bounds error";
   		WriteBlock(enc, i*segment_size, (i+1)*segment_size, v);
    }

    void initialize() {
    	parlay::parallel_for(0, aux_size, [&] (uint32_t i) {
    		write(i, aux[i]);
    	});
    }

    void restore() {
    	parlay::parallel_for(0, aux_size, [&] (uint32_t i) {
    		aux[i] = read(i);
    	});
    	parlay::parallel_for(0, enc.size()/2, [&] (uint32_t i) {
    		uint32_t idx1 = 2*i;
            uint32_t idx2 = 2*i + 1;
            if (enc[idx1] > enc[idx2]) {
            std::swap(enc[idx1], enc[idx2]);
            }
    	});
    }
};
