#ifndef BSP_HEADER
#define BSP_HEADER

#include "util/types.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

#include "xtensor/xarray.hpp"

namespace osiris {

///@brief simple generic container that paritions a hyber-box of arbitrary
/// dimensional space into nested binary sub-volumes, down to a fixed depth,
/// associating a T instance with each of the deepest sub-volumes. Splits along
/// dimensions in order they're provided by bounds_left and right, cycling back
/// the first dimension if/when the number of layers exceeds the number of
/// dimensions.
template <class T, class Point = xt::xarray<real>> class BinarySPTree {
public:
  /// @brief number of layers, number of partitions goes as 2^depth, partition
  /// volumes go with (1/2)^depth
  const int depth{};
  /// @brief number of dimensions of input space
  const int dimensions{};
  /// @brief outer bounds/corners of hyper-box in input space
  const Point bounds_left, bounds_right{};

  BinarySPTree(int depth, Point bounds_left, Point bounds_right,
               std::vector<T> data)
      : depth(depth--), dimensions(bounds_left.size()),
        bounds_left(bounds_left), bounds_right(bounds_right), data(data),
        root(std::make_unique<Node>(depth--, 0, bounds_left, bounds_right, 0,
                                    data.size() - 1)) {
    if (bounds_left.size() != bounds_right.size())
      throw std::runtime_error(
          "dimensions mismatch in bounds_left and bounds_right");
    for (int i = 0; i < dimensions; ++i)
      if (bounds_left[i] >= bounds_right[i]) {
        throw std::runtime_error(
            "bounds_left must be strictly less than bounds_right");
      }
  }

  /// @brief get reference to element associated with partition in which point
  /// resides
  T &operator[](const Point &p);
  /// @brief get const reference to element associated with partition in which
  /// point resides
  const T &at(const Point &p) const;

  using idx = int;

  auto begin() { return data.begin(); }
  auto end() { return data.end(); }
  auto cbegin() const { return data.cbegin(); }
  auto cend() const { return data.cend(); }

  int size() const { return std::pow(2, depth); }

private:
  /// @brief Encapsulates a partitioning of the input space into two
  /// sub-volumes, left and right
  class Node {
    std::unique_ptr<Node> left, right;

  protected:
    struct Split {
      std::pair<idx, idx> left;
      std::pair<idx, idx> right;
    };

    static Split split_down_the_middle(idx data_begin, idx data_end);

  public:
    ///@brief remaining depth to leaves, 0 indicating this node is a leaf P
    const int depth{};
    /// @brief dimension we split on
    const int split_dimension;
    /// @brief the boundary between the "left" and "right" partitions of
    /// spit_dimension
    const real bound;

    Node(int depth, int split_dimension, Point bounds_left, Point bounds_right,
         idx data_begin, idx data_end);

    /// @brief
    virtual idx operator[](const Point &point);
  };

  /// @brief
  class Leaf : public Node {
  private:
    const idx idx_left{}, idx_right{};

  public:
    Leaf(int next_split_dimension, Point bounds_left, Point bounds_right,
         BinarySPTree<T, Point>::idx idx_left,
         BinarySPTree<T, Point>::idx idx_right);

    idx operator[](const Point &point) final;
  };

  std::vector<T> data{};
  std::unique_ptr<Node> root{};
  bool is_within_bounds(const Point &point) const;
};

template <class T, class Point>
typename BinarySPTree<T, Point>::Node::Split
BinarySPTree<T, Point>::Node::split_down_the_middle(idx begin, idx end) {
  assert((end - begin + 1) % 2 == 0);
  auto size = (end - begin + 1) / 2;
  auto middle = begin + size - 1;
  return {{begin, middle}, {middle + 1, end}};
}

template <class T, class Point>
BinarySPTree<T, Point>::Node::Node(int depth, int split_dimension,
                                   Point bounds_left, Point bounds_right,
                                   idx idx_begin, idx idx_end)
    : depth(depth), split_dimension(split_dimension),
      bound((bounds_left[split_dimension] + bounds_right[split_dimension]) /
            2) {

  if (depth != 0) {
    auto size = idx_end - idx_begin + 1;
    auto expected_size_below = std::pow(2, depth + 1);
    assert(expected_size_below == size);
  } else {
    assert(idx_end - idx_begin + 1 == 2);
  }

  // determine next split dimension
  auto next_split_dimension = split_dimension + 1;
  if (next_split_dimension == bounds_left.size())
    next_split_dimension = 0;

  // construct new boundaries
  auto new_bounds_right = bounds_right;
  new_bounds_right[split_dimension] = bound;
  auto new_bounds_left = bounds_left;
  new_bounds_left[split_dimension] = bound;

  // determine remaining depth from daughter node to leaves
  const auto daughter_depths = depth - 1;

  if (depth > 1) {

    auto [left_data, right_data] = split_down_the_middle(idx_begin, idx_end);

    left = std::move(std::make_unique<BinarySPTree<T, Point>::Node>(
        daughter_depths, next_split_dimension, bounds_left, new_bounds_right,
        left_data.first, left_data.second));

    right = std::move(std::make_unique<BinarySPTree<T, Point>::Node>(
        daughter_depths, next_split_dimension, new_bounds_left, bounds_right,
        right_data.first, right_data.second));
  } else if (depth == 1) {
    assert(idx_end - idx_begin + 1 == 4);
    auto left_start = idx_begin;
    auto left_end = idx_begin + 1;
    auto right_start = idx_begin + 2;
    auto right_end = idx_begin + 3;

    left = std::move(std::make_unique<BinarySPTree<T, Point>::Leaf>(
        next_split_dimension, bounds_left, new_bounds_right, left_start,
        left_end));
    right = std::move(std::make_unique<BinarySPTree<T, Point>::Leaf>(
        next_split_dimension, new_bounds_left, bounds_right, right_start,
        right_end));
  }
}

template <class T, class Point>
typename BinarySPTree<T, Point>::idx
BinarySPTree<T, Point>::Node::operator[](const Point &point) {
  if (point(split_dimension) >= bound)
    return right->operator[](point);
  return left->operator[](point);
}

template <class T, class Point>
BinarySPTree<T, Point>::Leaf::Leaf(int split_dimension, Point bounds_left,
                                   Point bounds_right,
                                   BinarySPTree<T, Point>::idx idx_left,
                                   BinarySPTree<T, Point>::idx idx_right)
    : BinarySPTree<T, Point>::Node(0, split_dimension, bounds_left,
                                   bounds_right, idx_left, idx_right),
      idx_left(idx_left), idx_right(idx_right) {}

template <class T, class Point>
typename BinarySPTree<T, Point>::idx
BinarySPTree<T, Point>::Leaf::operator[](const Point &point) {
  if (point(this->split_dimension) >= this->bound)
    return idx_right;
  return idx_left;
}

template <class T, class Point>
T &BinarySPTree<T, Point>::operator[](const Point &point) {
  assert(point.size() == dimensions);
  if (not is_within_bounds(point)) {
    throw std::runtime_error("out of bounds");
  }
  return data[root->operator[](point)];
}

template <class T, class Point>
const T &BinarySPTree<T, Point>::at(const Point &point) const {
  assert(point.size() == dimensions);
  if (not is_within_bounds(point)) {
    throw std::runtime_error("out of bounds");
  }
  return data.at(root->at(point));
}

template <class T, class Point>
bool BinarySPTree<T, Point>::is_within_bounds(const Point &point) const {
  assert(point.size() == dimensions);
  for (int i = 0; i < dimensions; ++i) {
    if (point(i) < bounds_left(i) or point(i) > bounds_right(i))
      return false;
  }
  return true;
}

}; // namespace osiris

#endif
