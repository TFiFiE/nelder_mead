#ifndef NELDER_MEAD_HPP
#define NELDER_MEAD_HPP

#include <algorithm>
#include <array>
#include <cassert>
#include <random>
#include <vector>

template<class Number, unsigned int dimensions> struct NelderMead {
  typedef std::array<Number, dimensions> Input;
  typedef std::pair<Input, Number> Vertex;

  enum Step {
    INITIALIZATION,
    REFLECTION,
    EXPANSION,
    OUTSIDE_CONTRACTION,
    INSIDE_CONTRACTION,
    SHRINK
  };

  std::vector<Vertex> vertices;
  const typename std::vector<Vertex>::iterator best, pastWorst, secondBest, worst, secondWorst;
  const Number oneOverNumBest;
  Input bestCentroid;

  template<class Iterator>
  NelderMead(const Iterator begin, const Iterator end) :
    vertices{begin, end},
    best(std::begin(vertices)),
    pastWorst(std::end(vertices)),
    secondBest(next(best)),
    worst(prev(pastWorst)),
    secondWorst(prev(worst)),
    oneOverNumBest(Number(1) / (distance(best, worst))),
    bestCentroid(centroid(best, worst, oneOverNumBest))
  {
  }

  template<class Container> NelderMead(const Container& container) : NelderMead(begin(container), end(container)) {}

  static Input& input(Vertex& vertex) { return vertex.first; }

  static const Input& input(const Vertex& vertex) { return vertex.first; }

  static Number& output(Vertex& vertex) { return vertex.second; }

  static const Number& output(const Vertex& vertex) { return vertex.second; }

  static bool isBetterThan(const Vertex& lhs, const Vertex& rhs) { return output(lhs) < output(rhs); }

  template<class Iterator> static Input centroid(const Iterator begin, const Iterator end, const Number oneOverDistance)
  {
    Input result = {{0}};
    for (Iterator vertex = begin; vertex != end; ++vertex)
      for (unsigned int c = 0; c < dimensions; ++c)
        result[c] += input(*vertex)[c];
    for (unsigned int c = 0; c < dimensions; ++c)
      result[c] *= oneOverDistance;
    return result;
  }

  template<class Function> static Vertex extrapolate(const Input& origin, const Vertex& reference, const Number factor, Function function, const Step step)
  {
    Vertex result;
    for (unsigned int c = 0; c < dimensions; ++c)
      input(result)[c] = origin[c] + factor * (input(reference)[c] - origin[c]);
    output(result) = function(input(result), step);
    return result;
  }

  void updateBestCentroid(const Vertex& replacer)
  {
    for (unsigned int c = 0; c < dimensions; ++c)
      bestCentroid[c] += (input(replacer)[c] - input(*worst)[c]) * oneOverNumBest;
  }

  template<class Function>
  std::pair<Step, int> iteration(Function function,
                                 const Number inverseReflectionParameter = -1,
                                 const Number expansionParameter = 1 + 2 / Number(dimensions),
                                 const Number contractionParameter = 3 / Number(4) - 1 / Number(2 * dimensions),
                                 const Number shrinkageParameter = 1 - 1 / Number(dimensions))
  {
    assert(is_sorted(best, pastWorst, isBetterThan));

    const Vertex reflection = extrapolate(bestCentroid, *worst, inverseReflectionParameter, function, REFLECTION);

    if (isBetterThan(reflection, *secondWorst)) {
      if (isBetterThan(reflection, *best)) {
        const Vertex expansion = extrapolate(bestCentroid, reflection, expansionParameter, function, EXPANSION);

        move_backward(best, worst, pastWorst);
        if (isBetterThan(expansion, reflection)) {
          *best = expansion;
          updateBestCentroid(expansion);
          return {EXPANSION, 0};
        } else {
          *best = reflection;
          updateBestCentroid(reflection);
          return {REFLECTION, 0};
        }
      } else {
        const auto displaced = upper_bound(secondBest, secondWorst, reflection, isBetterThan);
        move_backward(displaced, worst, pastWorst);
        *displaced = reflection;
        updateBestCentroid(reflection);
        return {REFLECTION, distance(best, displaced)};
      }
    } else {
      if (isBetterThan(reflection, *worst)) {
        const Vertex outsideContraction = extrapolate(bestCentroid, reflection, contractionParameter, function, OUTSIDE_CONTRACTION);

        if (!isBetterThan(reflection, outsideContraction))
          return {OUTSIDE_CONTRACTION, contract(outsideContraction)};
      } else {
        const Vertex insideContraction = extrapolate(bestCentroid, *worst, contractionParameter, function, INSIDE_CONTRACTION);

        if (isBetterThan(insideContraction, *worst))
          return {INSIDE_CONTRACTION, contract(insideContraction)};
      }

      const Input oldBest = input(*best);
      for (auto vertex = secondBest; vertex != pastWorst; ++vertex)
        *vertex = extrapolate(input(*best), *vertex, shrinkageParameter, function, SHRINK);
      stable_sort(best, pastWorst, isBetterThan);
      bestCentroid = centroid(best, worst, oneOverNumBest);
      return {SHRINK, oldBest == input(*best) ? 1 : 0};
    }
  }

  int contract(const Vertex& contraction)
  {
    const auto displaced = upper_bound(best, worst, contraction, isBetterThan);
    move_backward(displaced, worst, pastWorst);
    *displaced = contraction;

    if (displaced != worst)
      updateBestCentroid(contraction);
    return distance(best, displaced);
  }

  template<class Function> void recalculate(Function function)
  {
    for (auto& vertex : vertices)
      output(vertex) = function(input(vertex), INITIALIZATION);
    stable_sort(best, pastWorst, isBetterThan);
    bestCentroid = centroid(best, worst, oneOverNumBest);
  }

  template<class Iterator, class Function> static NelderMead<Number, dimensions> create(const Iterator begin, const Iterator end, Function function)
  {
    const unsigned int numVertices = distance(begin, end);
    assert(numVertices > 1);
    std::vector<Vertex> vertices;
    vertices.reserve(numVertices);
    for (Iterator input = begin; input != end; ++input)
      vertices.emplace_back(*input, function(*input, INITIALIZATION));

    const auto simplexEnd = vertices.begin() + std::min(dimensions + 1, numVertices);
    partial_sort(vertices.begin(), simplexEnd, vertices.end(), isBetterThan);
    return NelderMead<Number, dimensions>(vertices.begin(), simplexEnd);
  }

  template<class Container, class Function> static NelderMead<Number, dimensions> create(const Container& container, Function function)
  {
    return create(begin(container), end(container), function);
  }

  template<class Generator>
  static std::vector<Input> randomPolytope(Generator& generator,
                                           const std::array<std::pair<Number, Number>, dimensions>& ranges,
                                           const unsigned int numVertices = dimensions + 1)
  {
    std::vector<Input> result(numVertices);
    for (unsigned int c = 0; c < dimensions; ++c) {
      std::uniform_real_distribution<> distribution(ranges[c].first, ranges[c].second);
      for (auto& vertex : result)
        vertex[c] = distribution(generator);
    }
    return result;
  }
};

#endif
