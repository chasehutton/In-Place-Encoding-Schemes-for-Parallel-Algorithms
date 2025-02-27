#include <iostream>
#include <cstdint>
#include <random>
#include <algorithm>
#include <set>
#include <chrono>
#include <vector>
#include <utility>
#include <numeric>
#include <functional>
#include <fstream>


#include "merge_random.h"
#include "merge_random2.h"
#include "block_size.h"


#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"

parlay::sequence<uint32_t> GenAlmostSorted(uint32_t n, uint32_t k) {
  assert(k <= n/2);
  auto R = parlay::tabulate(n, [&] (uint32_t i) {
    return i;
  });

  std::random_device rd;
  std::size_t seed = (static_cast<std::size_t>(rd()) << 32) ^ rd();
  parlay::random_generator gen(seed);
  std::uniform_int_distribution<uint64_t> dis(0, n/2 - 1);

  auto T = parlay::tabulate(k, [&] (uint32_t i) {
    auto r = gen[i];
    return std::make_tuple<uint32_t, uint32_t>(dis(r), dis(r));
  });

  for (auto&& t : T) {
    std::swap(R[std::get<0>(t)], R[std::get<1>(t) + n/2]);
  }

  parlay::integer_sort_inplace(parlay::make_slice(R.begin(), R.begin() + n/2));
  parlay::integer_sort_inplace(parlay::make_slice(R.begin() + n/2, R.end()));

  return R;
}

void ParReverse(parlay::slice<uint32_t *, uint32_t *> A) {
  auto Asize = A.size();
  if (Asize <= 2048) {
    for (int i = 0; i < std::floor(Asize / 2); i++) {
      std::swap(A[i], A[Asize - i - 1]);
    }
  } else {
    parlay::parallel_for(0, std::floor(Asize / 2), [&] (uint32_t i) {
      std::swap(A[i], A[Asize - i - 1]);
    });
  }
}

void ParReverseSplit(parlay::slice<uint32_t *, uint32_t *> A, parlay::slice<uint32_t *, uint32_t *> B) {
  auto Asize = A.size();
  auto Bsize = B.size();
  auto n = Asize + Bsize;
  if (n <= 2048) {
    for (int i = 0; i < std::floor((n) / 2); i++) {
      if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
      else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
      else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
    }
  } else {
    parlay::parallel_for(0, std::floor(n / 2), [&] (uint32_t i) {
      if (i < Asize && n - i - 1 < Asize) std::swap(A[i], A[n - i - 1]);
      else if (i < Asize && n - i - 1 >= Asize) std::swap(A[i], B[n - i - 1 - Asize]);
      else std::swap(B[i - Asize], B[n - i - 1 - Asize]);
    });
  }
}

// void StrongPIPMerge(parlay::sequence<uint32_t>& A, parlay::sequence<uint32_t>& B, uint32_t startA, uint32_t endA,
//                     uint32_t startB, uint32_t endB, uint32_t threshold) {
//   auto Asize = endA - startA;
//   auto Bsize = endB - startB;
//   if (Asize + Bsize <= threshold) {
//     merge(parlay::make_slice(A.begin() + startA, A.begin() + endA), parlay::make_slice(B.begin() + startB, B.begin() + endB));
//     return;
//   }

//   if (Asize == 0 || Bsize == 0) return;

//   if (Asize > Bsize) return StrongPIPMerge(B, A, startB, endB, startA, endA, threshold);

//   auto size = Asize + Bsize;
//   auto half = (size + 1) / 2;
//   auto left = startA;
//   auto right = endA;
//   uint32_t i;
//   uint32_t j;

//   while (left <= right) {
//     i = (left + right) / 2;
//     j = half - i + startB;

//     auto AleftMax  = (i == startA) ? std::numeric_limits<uint32_t>::min() : A[i - 1];
//     auto ArightMin = (i == endA) ? std::numeric_limits<uint32_t>::max() : A[i];

//     int BleftMax  = (j == startB) ? std::numeric_limits<uint32_t>::min() : B[j - 1];
//     int BrightMin = (j == endB) ? std::numeric_limits<uint32_t>::max() : B[j];

//     if (AleftMax <= BrightMin && BleftMax <= ArightMin) {
//       auto D1 = parlay::make_slice(A.begin() + startA + i, A.begin() + endA);
//       auto D2 = parlay::make_slice(B.begin() + startB, B.begin() + startB + j);
//       ParReverse(D1);
//       ParReverse(D2);
//       ParReverseSplit(D1,D2);
//       break;
//     } else if (AleftMax > BrightMin) {
//       right = i-1;
//     } else {
//       left = i+1;
//     }
//   }

//   parlay::par_do( 
//       [&] {
//       },
//       [&] {
          
//       }
//   );





// }

void GenRandom (uint64_t n, parlay::sequence<uint32_t>& seq, parlay::sequence<uint32_t>& seq32) {
    size_t mid = n / 2;

    parlay::parallel_do(
        [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin(), seq.begin() + mid)); },
        [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin() + mid, seq.end())); }
    );

    parlay::parallel_for(0, seq.size(), [&] (uint32_t i) {seq32[i] = static_cast<uint32_t>(seq[i]);});

    return;
}


// // This file contains GPT generated code for testing purposes

bool IsSorted(const parlay::sequence<uint32_t>& seq) {
  // A sequence with 0 or 1 elements is trivially sorted
  if (seq.size() <= 1) return true;

  // Check each adjacent pair
  for (size_t i = 1; i < seq.size(); i++) {
    if (seq[i] < seq[i - 1]) {
      return false;  // Found a descending pair
    }
  }
  return true;
}


bool IsSorted(const parlay::sequence<uint32_t>& seq, uint32_t s, uint32_t e) {

  // Check each adjacent pair
  for (size_t i = 1 + s; i < e; i++) {
    if (seq[i] < seq[i - 1]) {
      return false;  // Found a descending pair
    }
  }
  return true;
}


// void SaveSequenceToFile(parlay::sequence<uint32_t>& seq, const std::string& filename) {
//     std::ofstream outFile(filename);
//     if (!outFile) {
//         std::cerr << "Error: Unable to open file " << filename << " for writing.\n";
//         return;
//     }
    
//     outFile << "Failed Sequence:\n";
//     for (size_t i = 0; i < seq.size(); i++) {
//         outFile << seq[i] << (i % 10 == 9 ? "\n" : " ");
//     }
//     outFile.close();
//     std::cout << "Failed sequence written to " << filename << "\n";
// }



void driver(uint32_t n, uint32_t k) {
  parlay::sequence<uint32_t> times1(k);
  parlay::sequence<uint32_t> times2(k);
  parlay::sequence<uint32_t> times3(k);

	parlay::sequence<uint32_t> testSequence(n);
	std::random_device rd;
  uint64_t salt = rd();
	auto seq = parlay::random_permutation(n, salt);
  GenRandom(n, seq, testSequence);
    
    for (int i = 0; i < k; i++) {
      //testSequence = GenAlmostSorted(n, 1000);
      auto half = testSequence.size() / 2;
      auto A = testSequence.subseq(0, half);
      auto B = testSequence.subseq(half, testSequence.size());
      auto A1 = testSequence.subseq(0, half);
      auto B1 = testSequence.subseq(half, testSequence.size());
      auto A2 = testSequence.subseq(0, half);
      auto B2 = testSequence.subseq(half, testSequence.size());

      
      auto start3 = std::chrono::high_resolution_clock::now();
      version2::Merge2(A2, B2, b);
      auto end3 = std::chrono::high_resolution_clock::now();
      auto time3 = std::chrono::duration_cast<std::chrono::microseconds>(end3-start3);

      times3.push_back(time3.count());

      auto start2 = std::chrono::system_clock::now();
      auto R = parlay::merge(A1, B1);
      auto end2 = std::chrono::system_clock::now();
      auto time2 = std::chrono::duration_cast<std::chrono::microseconds>(end2-start2);

      // std::cout << "Time for Parlay Merge in iteration " << i << ": " << time2.count() << "\n";

      times2.push_back(time2.count());

      auto start1 = std::chrono::system_clock::now();
      Merge(A, B);
      auto end1 = std::chrono::system_clock::now();
      auto time1 = std::chrono::duration_cast<std::chrono::microseconds>(end1-start1);

      // std::cout << "Time for In-place Merge in iteration " << i << ": " << time1.count() << "\n";

      times1.push_back(time1.count());


    //   for (size_t i = 0; i < testSequence.size(); i++) {
    //     if (testSequence[i] != R[i]) {
    //       std::cout << "ERROR: Mismatch at index " << i 
    //                 << " => " << testSequence[i] << " vs. " << R[i] << "\n";
    //     }
    // }

      // if (!IsSorted(parlay::append(A,B))) {
      //   std::cout << "\nIteration " << i << " is not sorted! NEW \n";
      // }

      // if (!IsSorted(parlay::append(A2,B2))) {
      //   std::cout << "\nIteration " << i << " is not sorted! OLD\n";
      // }
    }

    auto t1 = parlay::reduce(times1) - times1[0];
    auto t2 = parlay::reduce(times2) - times2[0];
    auto t3 = parlay::reduce(times3) - times3[0];

    double avg_time1 = t1 / k;

    double avg_time2 = t2 / k;

    double avg_time3 = t3 / k;
    
    std::cout << "\nAvg Time In Microseconds for Old In-Place Merge: " << avg_time3 << "\n";
    std::cout << "\n\nAvg Time In Microseconds for In-Place Merge: " << avg_time1 << "\n";
    std::cout << "\nAvg Time In Microseconds for Parlay Merge: " << avg_time2 << "\n";
    std::cout << "\nSpeeddown: " << avg_time1/avg_time2 <<  "\n\n";
}

int main(int argc, char* argv[]) {
    uint32_t h = static_cast<uint32_t>(std::atoi(argv[1]));
    uint32_t size = h*(static_cast<uint32_t>(std::atoi(argv[2])));
    uint32_t tests = static_cast<uint32_t>(std::atoi(argv[3]));
    // auto A = parlay::make_slice(seq.begin(), seq.begin() + seq.size()/2);
    // auto B = parlay::make_slice(seq.begin() + seq.size()/2, seq.end());

   
    //auto R = parlay::merge(A,B);  
    // Merge(seq, 2048);
    // bool found = false;
    // parlay::sequence<uint32_t> testSequence;
    // while (!found) {
    //   testSequence = Gen2(size);
    //   if (testSequence[size/2 - 1] < 4194304 && testSequence[size/2 - 1]  > 4194300) {
    //     found = true;
    //   }
    // }
    // Merge(testSequence, 1024);
    driver(size, tests);


    //driver(size, 5);

//     parlay::sequence<uint32_t> testSequence = Gen2(size);
//     auto half = testSequence.size() / 2;
//     auto A = testSequence.subseq(0, half);
//     auto B = testSequence.subseq(half, testSequence.size());
// 
//     auto X = testSequence.subseq(0, half);
//         auto Y = testSequence.subseq(half, testSequence.size());

    // std::cout << "Testing...\n\n";

    // auto start = std::chrono::high_resolution_clock::now();

    // Merge(testSequence, 2048);

    // auto end = std::chrono::high_resolution_clock::now();

    // auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);

    // std::cout << "Time taken by IN-PLACE Merge: "
    //      << time.count() << " microseconds" << std::endl;

	  // buffer_merge(A, B);
    

    // start = std::chrono::high_resolution_clock::now();

    // auto r = parlay::merge(X, Y);

    // end = std::chrono::high_resolution_clock::now();

    // time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);

    // std::cout << "Time taken by Merge: "
    //      << time.count() << " microseconds" << std::endl;


    // if (testSequence.size() != r.size()) {
    //     std::cout << "ERROR: Different sizes!\n";
    //     std::cout << r.size() << "\n";
    //     std::cout << testSequence.size();
    //     return 1;
    // }

    // if (!IsSorted(testSequence)) {
    //     std::cout << "Not Sorted\n\n\n";
    //     // for (int i = 0; i < 1024; i++) {
    //     //     std::cout << testSequence[i] << "         ";
    //     //     if (i % 10 == 0) std::cout << "\n\n";
    //     // }
    //     std::cout << "\n\n\n\n";

    //     // auto e = ComputeAllInversions(testSequence);

    //     // for (auto r : e) {
    //     //     std::cout << "Inversion at indices " << r.first << " and " << r.second << "\n";
    //     //     std::cout << "Inverted Pair = " << "( " << testSequence[r.first] << ", " << testSequence[r.second] << " )\n"; 
    //     // }
    // }

	
    // for (size_t i = 0; i < A.size(); i++) {
    //         if (A[i] != r[i]) {
    //         std::cout << "ERROR: Mismatch at index " << i 
    //                     << " => " << A[i] << " vs. " << r[i] << "\n";
    //         return 1;
    //         }
    //     }
    //     for (size_t i = 0; i < B.size(); i++) {
    //                 if (B[i] != r[A.size()+i]) {
    //                 std::cout << "ERROR: Mismatch at index " << i 
    //                             << " => " << B[i] << " vs. " << r[A.size()+i] << "\n";
    //                 return 1;
    //                 }
    //             }
    
    // //     std::cout << "SUCCESS: Buffer Merge result matches parlay::merge.\n";
    
    // //Check if they are identical
    // for (size_t i = 0; i < testSequence.size(); i++) {
    //     if (testSequence[i] != r[i]) {
    //     std::cout << "ERROR: Mismatch at index " << i 
    //                 << " => " << testSequence[i] << " vs. " << r[i] << "\n";
    //     return 1;
    //     }
    // }

    // std::cout << "SUCCESS: In-place Merge result matches parlay::merge.\n";
    return 0;
}
