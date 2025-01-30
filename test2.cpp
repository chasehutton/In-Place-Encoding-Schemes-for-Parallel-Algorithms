#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

using namespace std;

// Function to find and print first inversion for each index
void findFirstInversions(const vector<int>& arr) {
    for (size_t i = 0; i < arr.size(); ++i) {
        for (size_t j = i + 1; j < arr.size(); ++j) {
            if (arr[i] > arr[j]) {
                cout << "First inversion for index " << i << ": (" << i << ", " << j << ") -> (" << arr[i] << ", " << arr[j] << ")\n";
                break; // Stop after finding the first inversion for index i
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>\n";
        return 1;
    }

    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Error opening file: " << argv[1] << "\n";
        return 1;
    }

    string line;
    bool firstLine = true;
    vector<int> numbers;
    int it = 0;
    while (getline(infile, line) && it < 200000) {
        if (firstLine) {
            firstLine = false; // Skip the header line
            continue;
        }
        
        stringstream ss(line);
        int num;
        while (ss >> num) {
            numbers.push_back(num);
        }

        it++;
    }
    
    if (!numbers.empty()) {
        //cout << "Processing entire sequence\n";
        findFirstInversions(numbers);
    }

    return 0;
}
