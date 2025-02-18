/*
 Case of vector
include "parallel_for_each.h"
#include <vector>
#include <iostream>

int main() {
  std::vector<int> data = {1, 2, 3, 4, 5};

  // Define a lambda function to apply to each element
  auto square = [](int& x) { x = x * x; };

  // Apply the function in parallel
  parallel::for_each(data.begin(), data.end(), square);

  // Print the results
  for (const auto& x : data) {
    std::cout << x << " ";
  }
  std::cout << std::endl;

  return 0;
}
*/

//Case to handle doubles

#include "parallel_for_each.h"
#include <iostream>

void example() {
    double arr[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

    parallel::for_each(arr, arr + 5, [] __host__ __device__ (double& x) {
        x *= 2.0;
    });

    for (double x : arr) {
        std::cout << x << " ";  // Should print: 2 4 6 8 10
    }
    std::cout << std::endl;
}

int main() {
    example();
    return 0;
}
