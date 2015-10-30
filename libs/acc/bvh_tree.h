/*
 * Copyright (C) 2015, Nils Moehrle
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#ifndef ACC_BVHTREE_HEADER
#define ACC_BVHTREE_HEADER

#include <deque>
#include <stack>
#include <algorithm>

#include <math/vector.h>

#include "primitives.h"

ACC_NAMESPACE_BEGIN

class BVHTree {
public:
    struct Hit {
        /* Parameter of the ray (distance of hit location). */
        float t;
        /* Index of the struck triangle. */
        std::size_t idx;
        /* Barycentric coordinates of hit location w.r.t. the triangle. */
        math::Vec3f bcoords;
    };

private:
    struct Node {
        std::size_t first;
        std::size_t last;
        Node* left;
        Node* right;
        AABB aabb;
        Node(std::size_t first, std::size_t last)
            : first(first), last(last), left(nullptr), right(nullptr) {}
        bool is_leaf() const {return left == nullptr && right == nullptr;}
    };
    std::vector<std::size_t> indices;
    std::vector<Tri> tris;

    Node *root;

    bool intersect(Ray const & ray, Node const * node, Hit * hit) const;

public:
    ~BVHTree();

    /* Constructs the BVH tree using the Surface Area Heuristic as
     * published in
     * "On fast Construction of SAH-based Bounding Volume Hierarchies"
     * by Ingo Wald (IEEE Symposium on Interactive Ray Tracing 2007)
     *
     * The mesh should be given as triangle index list and
     * a vector containing the 3D positions. */
    BVHTree(std::vector<std::size_t> const & faces,
        std::vector<math::Vec3f> const & vertices);

    bool intersect(Ray ray, Hit * hit_ptr) const;
};

BVHTree::~BVHTree() {
    std::deque<Node*> q;
    q.push_back(root);
    while (!q.empty()) {
        Node *node = q.back(); q.pop_back();
        if (!node->is_leaf()) {
            q.push_back(node->left);
            q.push_back(node->right);
        }
        delete node;
    }
};

BVHTree::BVHTree(std::vector<std::size_t> const & faces,
    std::vector<math::Vec3f> const & vertices) {
    
    std::size_t num_faces = faces.size() / 3;
    std::vector<AABB> aabbs(num_faces);
    std::vector<AABB> right_aabbs(num_faces);
    std::vector<Tri> ttris(num_faces);
    root = new Node(0, num_faces);
    for (std::size_t i = 0; i < aabbs.size(); ++i) {
        ttris[i].a = vertices[faces[i * 3 + 0]];
        ttris[i].b = vertices[faces[i * 3 + 1]];
        ttris[i].c = vertices[faces[i * 3 + 2]];

        calculate_aabb(ttris[i], &aabbs[i]);
        root->aabb += aabbs[i];
    }
    indices.resize(aabbs.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
        indices[i] = i;
    }

    std::deque<Node*> q;
    q.push_back(root);
    while (!q.empty()) {
        Node *node = q.back(); q.pop_back();
        std::size_t n = node->last - node->first;

        float min_cost = std::numeric_limits<float>::infinity();
        std::pair<std::size_t, std::size_t> split;
        for (std::size_t d = 0; d < 3; ++d) {
            std::sort(&indices[node->first], &indices[node->last],
                [&aabbs, d] (std::size_t first, std::size_t second) -> bool {
                    return mid(aabbs[first], d) < mid(aabbs[second], d)
                        || (mid(aabbs[first], d) == mid(aabbs[second], d)
                            && first < second);
                }
            );

            right_aabbs[node->last - 1] = aabbs[indices[node->last - 1]];
            for (std::size_t i = node->last - 1; i > node->first; --i) {
                right_aabbs[i - 1] = aabbs[indices[i - 1]] + right_aabbs[i];
            }
            node->aabb = right_aabbs[node->first];

            AABB left_aabb = aabbs[indices[node->first]];
            for (std::size_t i = node->first + 1; i < node->last; ++i) {
                std::size_t nl = i - node->first;
                std::size_t nr = n - nl;
                float cost = (surface_area(left_aabb) / surface_area(node->aabb) * nl
                + surface_area(right_aabbs[i]) / surface_area(node->aabb) * nr);
                if (cost <= min_cost) {
                    min_cost = cost;
                    split = std::make_pair(d, i);
                }

                left_aabb += aabbs[indices[i]];
            }
        }

        if (min_cost < n) {
            std::size_t d, i;
            std::tie(d, i) = split;
            if (d != 2) {
                std::sort(&indices[node->first], &indices[node->last],
                    [&aabbs, d] (std::size_t first, std::size_t second) -> bool {
                        return mid(aabbs[first], d) < mid(aabbs[second], d)
                            || (mid(aabbs[first], d) == mid(aabbs[second], d)
                                && first < second);
                    }
                );
            }

            node->left = new Node(node->first, i);
            node->right = new Node(i, node->last);
            q.push_back(node->left);
            q.push_back(node->right);
        }
    }

    tris.resize(ttris.size());
    for (std::size_t i = 0; i < indices.size(); ++i) {
        tris[i] = ttris[indices[i]];
    }
}

bool
BVHTree::intersect(Ray const & ray, Node const * node, Hit * hit) const {
    bool ret = false;
    for (std::size_t i = node->first; i < node->last; ++i) {
        float t;
        math::Vec3f bcoords;
        if (acc::intersect(ray, tris[i], &t, &bcoords)) {
            if (t > hit->t) continue;
            hit->idx = indices[i];
            hit->t = t;
            hit->bcoords = bcoords;
            ret = true;
        }
    }
    return ret;
}

bool
BVHTree::intersect(Ray ray, Hit * hit_ptr) const {
    Hit hit;
    hit.t = std::numeric_limits<float>::infinity();
    std::stack<Node const *> s;

    s.push(root);
    while (!s.empty()) {
        Node const *node = s.top(); s.pop();
        if (!node->is_leaf()) {
            float tmin_left, tmin_right;
            bool left = acc::intersect(ray, node->left->aabb, &tmin_left);
            bool right = acc::intersect(ray, node->right->aabb, &tmin_right);
            if (left && right) {
                if (tmin_left < tmin_right) {
                    s.push(node->right);
                    s.push(node->left);
                } else {
                    s.push(node->left);
                    s.push(node->right);
                }
            } else {
                if (right) s.push(node->right);
                if (left) s.push(node->left);
            }
        } else {
            if (intersect(ray, node, &hit)) {
                ray.tmax = hit.t;
            }
        }
    }

    *hit_ptr = hit;

    return hit.t < std::numeric_limits<float>::infinity();
}

ACC_NAMESPACE_END

#endif /* ACC_BVHTREE_HEADER */
