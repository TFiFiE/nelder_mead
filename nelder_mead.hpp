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
  const typename std::vector<Vertex>::iterator best, past_worst, second_best, worst, second_worst;
  const Number one_over_num_best;
  Input best_centroid;

  template<class Iterator>
  NelderMead(const Iterator begin, const Iterator end) :
    vertices{begin, end},
    best(std::begin(vertices)),
    past_worst(std::end(vertices)),
    second_best(std::next(best)),
    worst(std::prev(past_worst)),
    second_worst(std::prev(worst)),
    one_over_num_best(Number(1) / std::distance(best, worst)),
    best_centroid(centroid(best, worst, one_over_num_best))
  {
  }

  template<class Container> NelderMead(const Container& container) : NelderMead(std::begin(container), std::end(container)) {}

  static Input& input(Vertex& vertex) { return vertex.first; }

  static const Input& input(const Vertex& vertex) { return vertex.first; }

  static Number& output(Vertex& vertex) { return vertex.second; }

  static const Number& output(const Vertex& vertex) { return vertex.second; }

  static bool is_better_than(const Vertex& lhs, const Vertex& rhs) { return output(lhs) < output(rhs); }

  template<class Iterator> static Input centroid(const Iterator begin, const Iterator end, const Number one_over_distance)
  {
    Input result = {{0}};
    for (Iterator vertex = begin; vertex != end; ++vertex)
      for (unsigned int c = 0; c < dimensions; ++c)
        result[c] += input(*vertex)[c];
    for (unsigned int c = 0; c < dimensions; ++c)
      result[c] *= one_over_distance;
    return result;
  }

  template<class Function> static Vertex extrapolate(const Input& origin, const Vertex& reference, const Number factor, Function& function, const Step step)
  {
    Vertex result;
    for (unsigned int c = 0; c < dimensions; ++c)
      input(result)[c] = origin[c] + factor * (input(reference)[c] - origin[c]);
    output(result) = function(input(result), step);
    return result;
  }

  void update_best_centroid(const Vertex& replacer)
  {
    for (unsigned int c = 0; c < dimensions; ++c)
      best_centroid[c] += (input(replacer)[c] - input(*worst)[c]) * one_over_num_best;
  }

  template<class Function>
  std::pair<Step, int> iteration(Function& function,
                                 const Number inverse_reflection_parameter = -1,
                                 const Number expansion_parameter = 1 + 2 / Number(dimensions),
                                 const Number contraction_parameter = 3 / Number(4) - 1 / Number(2 * dimensions),
                                 const Number shrinkage_parameter = 1 - 1 / Number(dimensions))
  {
    assert(std::is_sorted(best, past_worst, is_better_than));

    const Vertex reflection = extrapolate(best_centroid, *worst, inverse_reflection_parameter, function, REFLECTION);

    if (is_better_than(reflection, *second_worst)) {
      if (is_better_than(reflection, *best)) {
        const Vertex expansion = extrapolate(best_centroid, reflection, expansion_parameter, function, EXPANSION);

        std::move_backward(best, worst, past_worst);
        if (is_better_than(expansion, reflection)) {
          *best = expansion;
          update_best_centroid(expansion);
          return {EXPANSION, 0};
        } else {
          *best = reflection;
          update_best_centroid(reflection);
          return {REFLECTION, 0};
        }
      } else {
        const auto displaced = std::upper_bound(second_best, second_worst, reflection, is_better_than);
        std::move_backward(displaced, worst, past_worst);
        *displaced = reflection;
        update_best_centroid(reflection);
        return {REFLECTION, std::distance(best, displaced)};
      }
    } else {
      if (is_better_than(reflection, *worst)) {
        const Vertex outside_contraction = extrapolate(best_centroid, reflection, contraction_parameter, function, OUTSIDE_CONTRACTION);

        if (!is_better_than(reflection, outside_contraction))
          return {OUTSIDE_CONTRACTION, contract(outside_contraction)};
      } else {
        const Vertex inside_contraction = extrapolate(best_centroid, *worst, contraction_parameter, function, INSIDE_CONTRACTION);

        if (is_better_than(inside_contraction, *worst))
          return {INSIDE_CONTRACTION, contract(inside_contraction)};
      }

      const Input old_best = input(*best);
      for (auto vertex = second_best; vertex != past_worst; ++vertex)
        *vertex = extrapolate(input(*best), *vertex, shrinkage_parameter, function, SHRINK);
      std::stable_sort(best, past_worst, is_better_than);
      best_centroid = centroid(best, worst, one_over_num_best);
      return {SHRINK, old_best == input(*best) ? 1 : 0};
    }
  }

  int contract(const Vertex& contraction)
  {
    const auto displaced = std::upper_bound(best, worst, contraction, is_better_than);
    std::move_backward(displaced, worst, past_worst);
    *displaced = contraction;

    if (displaced != worst)
      update_best_centroid(contraction);
    return std::distance(best, displaced);
  }

  template<class Function> void recalculate(Function& function)
  {
    for (auto& vertex : vertices)
      output(vertex) = function(input(vertex), INITIALIZATION);
    std::stable_sort(best, past_worst, is_better_than);
    best_centroid = centroid(best, worst, one_over_num_best);
  }

  template<class Iterator, class Function> static NelderMead<Number, dimensions> create(const Iterator begin, const Iterator end, Function& function)
  {
    const unsigned int num_vertices = std::distance(begin, end);
    assert(num_vertices > 1);
    std::vector<Vertex> vertices;
    vertices.reserve(num_vertices);
    for (Iterator input = begin; input != end; ++input)
      vertices.emplace_back(*input, function(*input, INITIALIZATION));

    const auto simplex_end = vertices.begin() + std::min(dimensions + 1, num_vertices);
    std::partial_sort(vertices.begin(), simplex_end, vertices.end(), is_better_than);
    return NelderMead<Number, dimensions>(vertices.begin(), simplex_end);
  }

  template<class Container, class Function> static NelderMead<Number, dimensions> create(const Container& container, Function& function)
  {
    return create(std::begin(container), std::end(container), function);
  }

  template<class Generator>
  static std::vector<Input> random_polytope(Generator& generator,
                                            const std::array<std::pair<Number, Number>, dimensions>& ranges,
                                            const unsigned int num_vertices = dimensions + 1)
  {
    std::vector<Input> result(num_vertices);
    for (unsigned int c = 0; c < dimensions; ++c) {
      std::uniform_real_distribution<> distribution(ranges[c].first, ranges[c].second);
      for (auto& vertex : result)
        vertex[c] = distribution(generator);
    }
    return result;
  }
};

#endif
