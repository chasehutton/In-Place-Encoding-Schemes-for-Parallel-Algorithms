#include <iostream>
#include <cstdint>

#include "utils.h"
#include "merge.h"



#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"


int main() {
    parlay::parallel_for(0, 100000, [&](size_t i) {
        int A[1000000] = {};
    });

    return 0;
}