#ifndef NELDER_MEAD_HPP
#define NELDER_MEAD_HPP

#include <algorithm>
#include <cassert>
#include <random>
#include <vector>

template<class Number> struct NelderMead {
  typedef std::vector<Number> Input;
  typedef std::pair<Input, Number> Vertex;

  enum Step {
    INITIALIZATION,
    REFLECTION,
    EXPANSION,
    OUTSIDE_CONTRACTION,
    INSIDE_CONTRACTION,
    SHRINK
  };

  const unsigned int num_dim;
  std::vector<Vertex> vertices;
  const typename std::vector<Vertex>::iterator best, past_worst, second_best, worst, second_worst;
  const Number one_over_num_best;
  Input best_centroid;

  template<class Container> NelderMead(const Container& container) : NelderMead(std::begin(container), std::end(container)) {}

  template<class Container> NelderMead(const Container& container, const unsigned int num_dim) : NelderMead(std::begin(container), std::end(container), num_dim)
  {
  }

  template<class Iterator> NelderMead(const Iterator begin, const Iterator end) : NelderMead(begin, end, begin->size()) {}

  template<class Iterator>
  NelderMead(const Iterator begin, const Iterator end, const unsigned int num_dim_) :
    num_dim(num_dim_),
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

  static Input& input(Vertex& vertex) { return vertex.first; }

  static const Input& input(const Vertex& vertex) { return vertex.first; }

  static Number& output(Vertex& vertex) { return vertex.second; }

  static const Number& output(const Vertex& vertex) { return vertex.second; }

  static bool is_better_than(const Vertex& lhs, const Vertex& rhs) { return output(lhs) < output(rhs); }

  static bool nan_is_better_than(const Vertex& lhs, const Vertex& rhs)
  {
    if (std::isnan(output(lhs)))
      return false;
    else if (std::isnan(output(rhs)))
      return true;
    else
      return is_better_than(lhs, rhs);
  }

  template<class Iterator> Input centroid(const Iterator begin, const Iterator end, const Number one_over_distance)
  {
    Input result(num_dim);
    for (Iterator vertex = begin; vertex != end; ++vertex)
      for (unsigned int c = 0; c < num_dim; ++c)
        result[c] += input(*vertex)[c];
    for (unsigned int c = 0; c < num_dim; ++c)
      result[c] *= one_over_distance;
    return result;
  }

  template<class Function> Vertex extrapolate(const Input& origin, const Vertex& reference, const Number factor, Function& function, const Step step)
  {
    Input result;
    result.reserve(num_dim);
    for (unsigned int c = 0; c < num_dim; ++c)
      result.emplace_back(origin[c] + factor * (input(reference)[c] - origin[c]));
    return {result, function(result, step)};
  }

  void update_best_centroid(const Vertex& replacer)
  {
    for (unsigned int c = 0; c < num_dim; ++c)
      best_centroid[c] += (input(replacer)[c] - input(*worst)[c]) * one_over_num_best;
  }

  template<class Function> std::pair<Step, int> iteration(Function& function)
  {
    const Number inverse_reflection_parameter = -1;
    const Number expansion_parameter = 1 + 2 / Number(num_dim);
    const Number contraction_parameter = 3 / Number(4) - 1 / Number(2 * num_dim);
    const Number shrinkage_parameter = 1 - 1 / Number(num_dim);
    return iteration(function, inverse_reflection_parameter, expansion_parameter, contraction_parameter, shrinkage_parameter);
  }

  template<class Function>
  std::pair<Step, int> iteration(Function& function,
                                 const Number inverse_reflection_parameter,
                                 const Number expansion_parameter,
                                 const Number contraction_parameter,
                                 const Number shrinkage_parameter)
  {
    assert(std::is_sorted(best, past_worst, nan_is_better_than));

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
      std::stable_sort(best, past_worst, nan_is_better_than);
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
    std::stable_sort(best, past_worst, nan_is_better_than);
    best_centroid = centroid(best, worst, one_over_num_best);
  }

  template<class Container, class Function> static NelderMead<Number> create(Function& function, const Container& container)
  {
    return create(function, std::begin(container), std::end(container));
  }

  template<class Container, class Function> static NelderMead<Number> create(Function& function, const Container& container, const unsigned int num_dim)
  {
    return create(function, std::begin(container), std::end(container), num_dim);
  }

  template<class Iterator, class Function> static NelderMead<Number> create(Function& function, const Iterator begin, const Iterator end)
  {
    return create(function, begin, end, begin->size());
  }

  template<class Iterator, class Function>
  static NelderMead<Number> create(Function& function, const Iterator begin, const Iterator end, const unsigned int num_dim)
  {
    const unsigned int num_vertices = std::distance(begin, end);
    assert(num_vertices > 1);
    std::vector<Vertex> vertices;
    vertices.reserve(num_vertices);
    for (Iterator input = begin; input != end; ++input)
      vertices.emplace_back(*input, function(*input, INITIALIZATION));

    const auto simplex_end = vertices.begin() + std::min(num_dim + 1, num_vertices);
    std::partial_sort(vertices.begin(), simplex_end, vertices.end(), nan_is_better_than);
    return NelderMead<Number>(vertices.begin(), simplex_end, num_dim);
  }

  template<class Generator> static std::vector<Input> random_polytope(Generator& generator, const std::vector<std::pair<Number, Number>>& ranges)
  {
    return random_polytope(generator, ranges, ranges.size() + 1);
  }

  template<class Generator>
  static std::vector<Input> random_polytope(Generator& generator, const std::vector<std::pair<Number, Number>>& ranges, const unsigned int num_vertices)
  {
    const unsigned int num_dim = ranges.size();
    std::vector<Input> result(num_vertices, Input(num_dim));
    for (unsigned int c = 0; c < num_dim; ++c) {
      std::uniform_real_distribution<> distribution(ranges[c].first, ranges[c].second);
      for (auto& vertex : result)
        vertex[c] = distribution(generator);
    }
    return result;
  }
};

#endif
