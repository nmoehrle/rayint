/*
 * Copyright (C) 2015, Nils Moehrle
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef ACC_KDTREE_HEADER
#define ACC_KDTREE_HEADER

#include <queue>
#include <stack>
#include <limits>
#include <atomic>

#include <math/vector.h>

#include "defines.h"

ACC_NAMESPACE_BEGIN

template <uint16_t K, typename IdxType = unsigned>
class KDTree {
public:
    #define NAI std::numeric_limits<IdxType>::max()
private:

    std::vector<math::Vector<float, K> > const & vertices;
    struct Node {
        typedef IdxType ID;
        decltype(K) d;
        IdxType first;
        IdxType last;
        IdxType vertex_id;
        Node::ID left;
        Node::ID right;
    };

    std::atomic<IdxType> num_nodes;
    std::vector<Node> nodes;
    typename Node::ID create_node(decltype(K) d, IdxType first, IdxType last) {
        typename Node::ID node_id = num_nodes++;
        Node & node = nodes[node_id];
        node.first = first;
        node.last = last;
        node.left = NAI;
        node.right = NAI;
        node.vertex_id = NAI;
        node.d = d;
        return node_id;
    }

public:
    KDTree(std::vector<math::Vector<float, K> > const & vertices);

    std::pair<IdxType, float>
    find_nn(math::Vector<float, K> point,
        float max_dist = std::numeric_limits<float>::infinity()) const;

    std::vector<std::pair<IdxType, float> >
    find_nns(math::Vector<float, K> point, std::size_t n,
        float max_dist = std::numeric_limits<float>::infinity()) const;
};

template <uint16_t K, typename IdxType>
KDTree<K, IdxType>::KDTree(std::vector<math::Vector<float, K> > const & vertices)
    : vertices(vertices), num_nodes(0) {

    std::size_t num_vertices = vertices.size();
    nodes.resize(num_vertices);

    std::vector<IdxType> indices(num_vertices);
    for (std::size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    std::deque<typename Node::ID> queue;
    queue.push_back(create_node(0, 0, num_vertices));
    while (!queue.empty()) {
        typename Node::ID node_id = queue.front(); queue.pop_front();
        Node & node = nodes[node_id];
        decltype(K) d = node.d;
        std::sort(&indices[node.first], &indices[node.last],
            [&vertices, d] (IdxType a, IdxType b) -> bool {
                return vertices[a][d] < vertices[b][d];
            }
        );
        d = (d + 1) % K;
        IdxType mid = (node.last + node.first) / 2;
        node.vertex_id = indices[mid];
        if (mid - node.first > 0) {
            node.left = create_node(d, node.first, mid);
            queue.push_back(node.left);
        }
        if (node.last - (mid + 1) > 0) {
            node.right = create_node(d, mid + 1, node.last);
            queue.push_back(node.right);
        }
    }
}

template <uint16_t K, typename IdxType>
std::pair<IdxType, float>
KDTree<K, IdxType>::find_nn(math::Vector<float, K> point, float max_dist) const {
    return find_nns(point, 1, max_dist)[0];
}

template <uint16_t K, typename IdxType>
std::vector<std::pair<IdxType, float> >
KDTree<K, IdxType>::find_nns(math::Vector<float, K> vertex, std::size_t n, float max_dist) const {

    std::pair<IdxType, float> nn = std::make_pair(NAI, max_dist);
    std::vector<std::pair<IdxType, float> > nns(n, nn);

    std::stack<std::pair<typename Node::ID, bool> > s;
    s.emplace(0, true);
    while (!s.empty()) {
        typename Node::ID node_id;
        bool down;
        std::tie(node_id, down) = s.top();
        s.pop();

        if (node_id == NAI) continue;

        Node const & node = nodes[node_id];

        float diff = vertex[node.d] - vertices[node.vertex_id][node.d];
        if (down) {
            float dist = (vertex - vertices[node.vertex_id]).norm();
            if (dist < max_dist) {
                nns.emplace_back(node.vertex_id, dist);
                std::sort(nns.begin(), nns.end(),
                    [] (std::pair<IdxType, float> a, std::pair<IdxType, float> b) -> bool {
                        return a.second < b.second;
                    }
                );
                nns.pop_back();
                max_dist = nns.back().second;
            }

            if (node.left == NAI && node.right == NAI) continue;

            s.emplace(node_id, false);
            if (diff < 0.0f) {
                s.emplace(node.left, true);
            } else {
                s.emplace(node.right, true);
            }
        } else {
            if (std::abs(diff) >= max_dist) continue;

            if (diff < 0.0f) {
                s.emplace(node.right, true);
            } else {
                s.emplace(node.left, true);
            }
        }
    }
    return nns;
}

ACC_NAMESPACE_END

#endif /* ACC_KDTREE_HEADER */
