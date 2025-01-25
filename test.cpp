#include <iostream>
#include <cstdint>
#include <random>
#include <algorithm>
#include <set>
#include <chrono>
#include <vector>
#include <utility>

#include "utils.h"
#include "merge.h"
#include "buffer.h"



#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include "parlay/primitives.h"

// This file contains GPT generated code for testing purposes



std::vector<std::pair<size_t,size_t>>
ComputeAllInversions(const parlay::sequence<uint32_t>& seq)
{
  std::vector<std::pair<size_t,size_t>> inversions;
  size_t n = seq.size();

  for (size_t i = 0; i < n; i++) {
    for (size_t j = i + 1; j < n; j++) {
      if (seq[i] > seq[j]) {
        // This is an inversion
        inversions.push_back({i, j});
      }
    }
  }
  return inversions;
}


// A helper function that merges two sorted halves [start,mid) and [mid,end)
// in-place, counting how many cross-inversions occur.
// Returns the number of inversions found in the merge step.
static std::size_t mergeAndCount(parlay::sequence<uint32_t>& arr,
                                 std::size_t start,
                                 std::size_t mid,
                                 std::size_t end)
{
  // We'll copy the two halves, then merge back
  std::size_t lenA = mid - start;
  std::size_t lenB = end - mid;

  parlay::sequence<uint32_t> left(lenA);
  parlay::sequence<uint32_t> right(lenB);

  // Copy data
  for (std::size_t i = 0; i < lenA; i++) {
    left[i] = arr[start + i];
  }
  for (std::size_t j = 0; j < lenB; j++) {
    right[j] = arr[mid + j];
  }

  // Merge them back into arr[start..end)
  std::size_t i = 0;
  std::size_t j = 0;
  std::size_t k = start;

  std::size_t inversions = 0;

  while (i < lenA && j < lenB) {
    if (left[i] <= right[j]) {
      arr[k++] = left[i++];
    } else {
      // left[i] > right[j] => each element left[i..(lenA-1)]
      // is greater than right[j]. So we add (lenA - i) to inversions.
      arr[k++] = right[j++];
      inversions += (lenA - i);
    }
  }

  // Copy any remaining
  while (i < lenA) {
    arr[k++] = left[i++];
  }
  while (j < lenB) {
    arr[k++] = right[j++];
  }

  return inversions;
}

// Recursive function that sorts arr[start..end) and returns #inversions
static std::size_t mergeSortCount(parlay::sequence<uint32_t>& arr,
                                  std::size_t start,
                                  std::size_t end)
{
  std::size_t n = end - start;
  if (n <= 1) {
    // 0 or 1 element => no inversions
    return 0;
  }
  // Divide
  std::size_t mid = start + (n / 2);

  // Conquer
  std::size_t leftCount = mergeSortCount(arr, start, mid);
  std::size_t rightCount = mergeSortCount(arr, mid, end);

  // Combine
  std::size_t mergeCount = mergeAndCount(arr, start, mid, end);

  return leftCount + rightCount + mergeCount;
}

// Public function: returns the total #inversions in the sequence
std::size_t CountInversions(parlay::sequence<uint32_t> seq) {
  // We'll sort a copy of seq in-place (or you could sort seq directly
  // if you don't mind altering it).
  return mergeSortCount(seq, 0, seq.size());
}

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


parlay::sequence<uint32_t> generateUniqueTwoSortedHalves(uint32_t size) {
    if (size % 2 != 0) {
        throw std::invalid_argument("Size must be even.");
    }
    size_t half = size / 2;

    // Random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<uint32_t> dis(1, 100000000);

    std::set<uint32_t> globalSet;
    while (globalSet.size() < size) {
        globalSet.insert(dis(gen));
    }

    parlay::sequence<uint32_t> allVals(globalSet.begin(), globalSet.end());
    std::shuffle(allVals.begin(), allVals.end(), gen);

    parlay::sequence<uint32_t> result(size);
    std::copy(allVals.begin(), allVals.begin() + half, result.begin());
    std::copy(allVals.begin() + half, allVals.end(),result.begin() + half);

    std::sort(result.begin(),  result.begin() + half);
    std::sort(result.begin() + half,  result.end());

    return result;
}

bool CheckEachHalfSorted(const parlay::sequence<uint32_t>& testSequence) {
  // Ensure testSequence has even length
  if (testSequence.size() % 2 != 0) {
    std::cerr << "CheckEachHalfSorted: sequence size is not even!\n";
    return false;
  }

  size_t n = testSequence.size();
  size_t half = n / 2;

  // 1) Check the first half: [0..(half-1)]
  for (size_t i = 1; i < half; i++) {
    if (testSequence[i] < testSequence[i - 1]) {
      std::cerr << "First half NOT sorted at index " << i
                << " (value " << testSequence[i]
                << " < " << testSequence[i - 1] << " at index " << (i - 1) << ")\n";
      return false;
    }
  }

  // 2) Check the second half: [half..(n-1)]
  for (size_t i = half + 1; i < n; i++) {
    if (testSequence[i] < testSequence[i - 1]) {
      std::cerr << "Second half NOT sorted at index " << i
                << " (value " << testSequence[i]
                << " < " << testSequence[i - 1] << " at index " << (i - 1) << ")\n";
      return false;
    }
  }

  // If we reach here, both halves are internally sorted
  return true;
}

int main() {
    size_t size = 4194304;

    parlay::sequence<uint32_t> testSequence = generateUniqueTwoSortedHalves(size);
    assert(CheckEachHalfSorted(testSequence));

    std::cout << "Testing...\n\n";

    auto start = std::chrono::high_resolution_clock::now();

    Merge(testSequence, 2048);

    auto end = std::chrono::high_resolution_clock::now();

    auto time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);

    std::cout << "Time taken by IN-PLACE Merge: "
         << time.count() << " microseconds" << std::endl;


    auto half = testSequence.size() / 2;
    auto A = testSequence.subseq(0, half);
    auto B = testSequence.subseq(half, testSequence.size());

    start = std::chrono::high_resolution_clock::now();

    auto r = parlay::merge(A, B);

    end = std::chrono::high_resolution_clock::now();

    time = std::chrono::duration_cast<std::chrono::microseconds>(end-start);

    std::cout << "Time taken by Merge: "
         << time.count() << " microseconds" << std::endl;


    if (testSequence.size() != r.size()) {
        std::cout << "ERROR: Different sizes!\n";
        std::cout << r.size() << "\n";
        std::cout << testSequence.size();
        return 1;
    }

    if (!IsSorted(testSequence)) {
        std::cout << "Not Sorted\n\n\n";
        std::cout << "Number of Inversion: " << CountInversions(testSequence);
        // for (int i = 0; i < 1024; i++) {
        //     std::cout << testSequence[i] << "         ";
        //     if (i % 10 == 0) std::cout << "\n\n";
        // }
        std::cout << "\n\n\n\n";

        // auto e = ComputeAllInversions(testSequence);

        // for (auto r : e) {
        //     std::cout << "Inversion at indices " << r.first << " and " << r.second << "\n";
        //     std::cout << "Inverted Pair = " << "( " << testSequence[r.first] << ", " << testSequence[r.second] << " )\n"; 
        // }
    }
    
    // Check if they are identical
    for (size_t i = 0; i < testSequence.size(); i++) {
        if (testSequence[i] != r[i]) {
        std::cout << "ERROR: Mismatch at index " << i 
                    << " => " << testSequence[i] << " vs. " << r[i] << "\n";
        return 1;
        }
    }

    std::cout << "SUCCESS: In-place Merge result matches parlay::merge.\n";
    return 0;
    return 0;
}
