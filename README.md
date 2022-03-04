# Example usage:
```C++
#include "nelder_mead.hpp"
#include <iostream>

template<class Number> inline Number square(const Number number) { return number * number; }

template<class Number> struct Rosenbrock {
  Number operator()(const std::array<Number, 2>& input,
                    const typename NelderMead<Number, 2>::Step) const
  {
    return square(1 - input[0]) + 100 * square(input[1] - square(input[0]));
  }
};

template<class Number> struct Himmelblau {
  Number operator()(const std::array<Number, 2>& input,
                    const typename NelderMead<Number, 2>::Step) const
  {
    return square(square(input[0]) + input[1] - 11) + square(input[0] + square(input[1]) - 7);
  }
};

template<class Number, class Generator, class Function>
void find_optimum(Generator& generator, Function function)
{
  typedef NelderMead<Number, 2> NelderMead;

  const auto simplex = NelderMead::random_polytope(generator, {{{-4, 4}, {-4, 4}}});
  auto nelder_mead = NelderMead::create(simplex, function);

  while (true) {
    const auto old_worst = *nelder_mead.worst;
    nelder_mead.iteration(function);
    if (old_worst == *nelder_mead.worst)
      break;
  }

  const auto result = *nelder_mead.best;
  std::cout << result.first[0] << ',' << result.first[1] << ':' << result.second << '\n';
}

int main()
{
  std::mt19937 generator((std::random_device())());
  typedef long double Number;
  find_optimum<Number>(generator, Rosenbrock<Number>());
  find_optimum<Number>(generator, Himmelblau<Number>());
}
```
