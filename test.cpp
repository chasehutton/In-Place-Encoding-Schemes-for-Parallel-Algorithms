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


#include "parlay/internal/merge.h"
#include "IPmerge.h"
#include "block_size.h"


#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"

void Gen (uint64_t n, parlay::sequence<uint32_t>& seq, parlay::sequence<uint32_t>& seq32) {
    size_t mid = n / 2;

    parlay::parallel_do(
        [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin(), seq.begin() + mid)); },
        [&] { parlay::integer_sort_inplace(parlay::make_slice(seq.begin() + mid, seq.end())); }
    );

    parlay::parallel_for(0, seq.size(), [&] (uint32_t i) {seq32[i] = static_cast<uint32_t>(seq[i]);});

    return;
}


bool IsSorted(const parlay::sequence<uint32_t>& seq) {
  for (size_t i = 1; i < seq.size(); i++) {
    if (seq[i] < seq[i - 1]) {
      return false; 
    }
  }
  return true;
}

void driver(uint32_t n, uint32_t k) {
  assert(k > 1);
  parlay::sequence<uint32_t> times1(k);
  parlay::sequence<uint32_t> times2(k);

	parlay::sequence<uint32_t> testSequence(n);
	std::random_device rd;
  uint64_t salt = rd();
	auto seq = parlay::random_permutation(n, salt);
  Gen(n, seq, testSequence);
    
    for (int i = 0; i < k; i++) {
      auto half = testSequence.size() / 2;
      auto A = testSequence.subseq(0, half);
      auto B = testSequence.subseq(half, testSequence.size());
      auto A1 = testSequence.subseq(0, half);
      auto B1 = testSequence.subseq(half, testSequence.size());

      auto start1 = std::chrono::system_clock::now();
      IPmerge(A, B);
      auto end1 = std::chrono::system_clock::now();
      auto time1 = std::chrono::duration_cast<std::chrono::microseconds>(end1-start1);
      times1.push_back(time1.count());

      auto start2 = std::chrono::system_clock::now();
      auto R = parlay::merge(A1, B1);
      auto end2 = std::chrono::system_clock::now();
      auto time2 = std::chrono::duration_cast<std::chrono::microseconds>(end2-start2);
      times2.push_back(time2.count());
    }

    auto t1 = parlay::reduce(times1) - times1[0];
    auto t2 = parlay::reduce(times2) - times2[0];

    double avg_time1 = t1 / (k-1);
    double avg_time2 = t2 / (k-1);
    
    std::cout << "\n\nAvg Time In Microseconds for In-Place Merge: " << avg_time1 << "\n";
    std::cout << "\nAvg Time In Microseconds for Parlay Merge: " << avg_time2 << "\n";
    std::cout << "\nSpeeddown: " << avg_time1/avg_time2 <<  "\n\n";
}

void generateSequencesAndWriteToFile(size_t n, const std::string& filename) {
  if (n % 2 != 0) {
    throw std::runtime_error("Error: n must be even.");
  }

  auto indices = parlay::tabulate(n, [&](size_t i) {
    return static_cast<uint32_t>(i);
  });


  std::random_device rd;
  uint64_t salt = rd();
  parlay::random rng(salt);
  parlay::random_shuffle(indices, rng);

  size_t half = n / 2;
  parlay::sequence<uint32_t> A(half), B(half);
  for (size_t i = 0; i < half; i++) {
    A[i] = indices[i];
  }
  for (size_t i = 0; i < half; i++) {
    B[i] = indices[half + i];
  }

  parlay::sort_inplace(A);
  parlay::sort_inplace(B);

  std::ofstream out(filename);
  if (!out.is_open()) {
    throw std::runtime_error("Error opening file for writing: " + filename);
  }

  out << half << "\n";  
  for (size_t i = 0; i < half; i++) {
    out << A[i] << "\n";
  }
  for (size_t i = 0; i < half; i++) {
    out << B[i] << "\n";
  }
  out.close();
}

void readSequencesAndMerge(const std::string& filename) {
  std::ifstream in(filename);
  if (!in.is_open()) {
    throw std::runtime_error("Error opening file for reading: " + filename);
  }

  size_t half;
  in >> half;

  parlay::sequence<uint32_t> A(half), B(half);
  for (size_t i = 0; i < half; i++) {
    in >> A[i];
  }
  for (size_t i = 0; i < half; i++) {
    in >> B[i];
  }
  in.close();

  auto C = parlay::merge(A, B);
  // do something
}


int main(int argc, char* argv[]) {
    uint32_t size =  b*(static_cast<uint32_t>(std::atoi(argv[1])));
    uint32_t tests = static_cast<uint32_t>(std::atoi(argv[2]));

    driver(size, tests);
    return 0;
}
