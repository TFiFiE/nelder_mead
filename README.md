# Example usage:
```
#include <iostream>
#include "nelder_mead.hpp"

template<class Number> inline Number square(const Number number) {return number*number;}

template<class Number>
struct Rosenbrock {
  Number operator()(const std::array<Number,2>& input,const typename NelderMead<Number,2>::Step) const
  {
    return square(1-input[0])+100*square(input[1]-square(input[0]));
  }
};

template<class Number>
struct Himmelblau {
  Number operator()(const std::array<Number,2>& input,const typename NelderMead<Number,2>::Step) const
  {
    return square(square(input[0])+input[1]-11)+square(input[0]+square(input[1])-7);
  }
};

template<class Number,class Generator,class Function>
void findOptimum(Generator& generator,Function function)
{
  typedef NelderMead<Number,2> NelderMead;

  const auto simplex=NelderMead::randomPolytope(generator,{{{-4,4},{-4,4}}});
  auto nelderMead=NelderMead::create(simplex,function);

  while (true) {
    const auto oldWorst=*nelderMead.worst;
    nelderMead.iteration(function);
    if (oldWorst==*nelderMead.worst)
      break;
  }

  const auto result=*nelderMead.best;
  std::cout<<result.first[0]<<','<<result.first[1]<<':'<<result.second<<'\n';
}

int main()
{
  std::mt19937 generator((std::random_device())());
  typedef long double Number;
  findOptimum<Number>(generator,Rosenbrock<Number>());
  findOptimum<Number>(generator,Himmelblau<Number>());
}
```
