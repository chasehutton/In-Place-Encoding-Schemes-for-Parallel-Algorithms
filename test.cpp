#include <iostream>
#include <cstdint>

#include "utils.h"
#include "merge.h"



#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"

void PrintSlice(parlay::slice<uint32_t*, uint32_t*> A) {
    for (int i = 0; i < A.size(); i++) {
        std::cout << A[i] << std::endl;
    }
}

int main() {
    uint32_t x = 100;
    auto seq = parlay::tabulate<uint32_t>(x, [&] (auto i) {
        return i;
    });

    PrintSlice(parlay::make_slice(seq.begin() + 5, seq.end() - 20));
    return 0;
}