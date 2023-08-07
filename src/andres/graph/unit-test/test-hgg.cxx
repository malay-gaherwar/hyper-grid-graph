#include <iostream>
#include <vector>
#include <array>

// Define the OffsetVector type alias with a template parameter
template <unsigned char D>
using OffsetVector = std::vector<std::array<size_t, D>>;

// Function to input and store the offset for a D-dimensional grid graph
template <unsigned char D>
OffsetVector<D> inputOffset() {
    OffsetVector<D> offset;
    size_t noffset;
    std::cout << "Enter the number of offset values for the D-dimensional grid graph:" << std::endl;
    std::cin >> noffset;

    for (size_t i = 0; i < noffset; ++i) {
        std::array<size_t, D> newArray;

        // Add elements to the array
        std::cout << "Enter the " << (i + 1) << " Offset list." << std::endl;
        for (size_t j = 0; j < D; ++j) {
            std::cout << "Enter the " << (j + 1) << " value:" << std::endl;
            size_t value;
            std::cin >> value;
            newArray[j] = value;
        }
        offset.push_back(newArray);
    }

    return offset;
}

int main() {
    const unsigned char D = 2; // Specify the dimensionality of the grid graph

    // Input the offset values based on the dimension D
    OffsetVector<D> offset = inputOffset<D>();

    // Print the stored offset values
    std::cout << "Offset values for the D-dimensional grid graph:" << std::endl;
    for (size_t i = 0; i < D; ++i) {
        for (const auto& value : offset[i]) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
