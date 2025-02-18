#include "parallel_for_each.h"
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
