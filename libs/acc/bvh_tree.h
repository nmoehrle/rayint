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
#include <cassert>
#include <algorithm>
#include <atomic>
#include <thread>

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
            : first(first), last(last), left(nullptr), right(nullptr),
            aabb({math::Vec3f(inf), math::Vec3f(-inf)}) {}
        bool is_leaf() const {return left == nullptr && right == nullptr;}
    };
    std::vector<std::size_t> indices;
    std::vector<Tri> tris;

    Node *root;

    std::pair<Node *, Node *> sbsplit(Node * node, std::vector<AABB> const & aabbs);
    std::pair<Node *, Node *> bsplit(Node * node, std::vector<AABB> const & aabbs);
    std::pair<Node *, Node *> ssplit(Node * node, std::vector<AABB> const & aabbs);
    void split(Node *node, std::vector<AABB> const & aabbs,
        std::atomic<int> * num_threads);

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
        std::vector<math::Vec3f> const & vertices,
        int max_threads = 2 * std::thread::hardware_concurrency());

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

#define NUM_BINS 64
struct Bin {
    std::size_t n;
    AABB aabb;
};


void BVHTree::split(Node *node, std::vector<AABB> const & aabbs,
        std::atomic<int> * num_threads) {

    Node *left, *right;
    if ((*num_threads -= 1) >= 1) {
        std::tie(left, right) = sbsplit(node, aabbs);
        if (left && right) {
            std::thread other(&BVHTree::split, this, left, std::cref(aabbs), num_threads);
            split(right, aabbs, num_threads);
            other.join();
        }
    } else {
        std::deque<Node*> queue;
        queue.push_back(node);
        while (!queue.empty()) {
            Node *node = queue.back(); queue.pop_back();
            std::tie(left, right) = sbsplit(node, aabbs);
            if (left && right) {
                queue.push_back(left);
                queue.push_back(right);
            }
        }
    }
    *num_threads += 1;
}

std::pair<BVHTree::Node *, BVHTree::Node *>
BVHTree::sbsplit(Node * node, std::vector<AABB> const & aabbs) {
    std::size_t n = node->last - node->first;
    if (n > NUM_BINS) {
        return bsplit(node, aabbs);
    } else {
        return ssplit(node, aabbs);
    }
}

std::pair<BVHTree::Node *, BVHTree::Node *>
BVHTree::bsplit(Node * node, std::vector<AABB> const & aabbs) {
    std::size_t n = node->last - node->first;

    std::array<Bin, NUM_BINS> bins;
    std::array<AABB, NUM_BINS> right_aabbs;
    std::vector<unsigned char> bin(n);

    float min_cost = std::numeric_limits<float>::infinity();
    std::pair<std::size_t, unsigned char> split;
    for (std::size_t d = 0; d < 3; ++d) {
        float min = node->aabb.min[d];
        float max = node->aabb.max[d];
        for (Bin & bin : bins) {
            bin = {0, {math::Vec3f(inf), math::Vec3f(-inf)}};
        }
        for (std::size_t i = node->first; i < node->last; ++i) {
            AABB const & aabb = aabbs[indices[i]];
            unsigned char idx = ((mid(aabb, d) - min) / (max - min)) * (NUM_BINS - 1);
            bins[idx].aabb += aabb;
            bins[idx].n += 1;
            bin[i - node->first] = idx;
        }

        right_aabbs[NUM_BINS - 1] = bins[NUM_BINS - 1].aabb;
        for (std::size_t i = NUM_BINS - 1; i > 0; --i) {
            right_aabbs[i - 1] = bins[i - 1].aabb + right_aabbs[i];
        }

        AABB left_aabb = bins[0].aabb;
        std::size_t nl = bins[0].n;
        for (std::size_t idx = 1; idx < NUM_BINS; ++idx) {
            std::size_t nr = n - nl;
            float cost = (surface_area(left_aabb) / surface_area(node->aabb) * nl
            + surface_area(right_aabbs[idx]) / surface_area(node->aabb) * nr);
            if (cost <= min_cost) {
                min_cost = cost;
                split = std::make_pair(d, idx);
            }

            nl += bins[idx].n;
            left_aabb += bins[idx].aabb;
        }
    }

    if (min_cost >= n) return std::make_pair(nullptr, nullptr);

    std::size_t d;
    unsigned char sidx;
    std::tie(d, sidx) = split;

    float min = node->aabb.min[d];
    float max = node->aabb.max[d];
    for (Bin & bin : bins) {
        bin = {0, {math::Vec3f(inf), math::Vec3f(-inf)}};
    }
    for (std::size_t i = node->first; i < node->last; ++i) {
        AABB const & aabb = aabbs[indices[i]];
        unsigned char idx = ((mid(aabb, d) - min) / (max - min)) * (NUM_BINS - 1);
        bins[idx].aabb += aabb;
        bins[idx].n += 1;
        bin[i - node->first] = idx;
    }

    std::size_t l = node->first;
    std::size_t r = node->last - 1;
    while (l < r) {
        if (bin[l - node->first] < sidx) {
            l += 1;
            continue;
        }
        if (bin[r - node->first] >= sidx) {
            r -= 1;
            continue;
        }
        std::swap(bin[l - node->first], bin[r - node->first]);
        std::swap(indices[l], indices[r]);
    }
    assert(l == r);
    std::size_t m = bin[(l&r) - node->first] >= sidx ? (l&r) : (l&r) + 1;

    node->left = new Node(node->first, m);
    node->right = new Node(m, node->last);
    for (std::size_t idx = 0; idx < NUM_BINS; ++idx) {
        if (idx < sidx) {
            node->left->aabb += bins[idx].aabb;
        } else {
            node->right->aabb += bins[idx].aabb;
        }
    }

    return std::make_pair(node->left, node->right);
}

std::pair<BVHTree::Node *, BVHTree::Node *>
BVHTree::ssplit(Node * node, std::vector<AABB> const & aabbs) {
    std::size_t n = node->last - node->first;

    float min_cost = std::numeric_limits<float>::infinity();
    std::pair<std::size_t, std::size_t> split;
    std::vector<AABB> right_aabbs(n);
    for (std::size_t d = 0; d < 3; ++d) {
        std::sort(&indices[node->first], &indices[node->last],
            [&aabbs, d] (std::size_t first, std::size_t second) -> bool {
                return mid(aabbs[first], d) < mid(aabbs[second], d)
                    || (mid(aabbs[first], d) == mid(aabbs[second], d)
                        && first < second);
            }
        );

        right_aabbs[n - 1] = aabbs[indices[node->last - 1]];
        for (std::size_t i = node->last - 1; i > node->first; --i) {
            right_aabbs[i - 1 - node->first] = aabbs[indices[i - 1]]
                + right_aabbs[i - node->first];
        }
        node->aabb = right_aabbs[0];

        AABB left_aabb = aabbs[indices[node->first]];
        for (std::size_t i = node->first + 1; i < node->last; ++i) {
            std::size_t nl = i - node->first;
            std::size_t nr = n - nl;
            float cost = (surface_area(left_aabb) / surface_area(node->aabb) * nl
            + surface_area(right_aabbs[nl]) / surface_area(node->aabb) * nr);
            if (cost <= min_cost) {
                min_cost = cost;
                split = std::make_pair(d, i);
            }

            left_aabb += aabbs[indices[i]];
        }
    }

    if (min_cost >= n) return std::make_pair(nullptr, nullptr);

    std::size_t d, i;
    std::tie(d, i) = split;
    std::sort(&indices[node->first], &indices[node->last],
        [&aabbs, d] (std::size_t first, std::size_t second) -> bool {
            return mid(aabbs[first], d) < mid(aabbs[second], d)
                || (mid(aabbs[first], d) == mid(aabbs[second], d)
                    && first < second);
        }
    );

    node->left = new Node(node->first, i);
    node->right = new Node(i, node->last);
    return std::make_pair(node->left, node->right);
}

BVHTree::BVHTree(std::vector<std::size_t> const & faces,
    std::vector<math::Vec3f> const & vertices, int max_threads) {

    std::size_t num_faces = faces.size() / 3;
    std::vector<AABB> aabbs(num_faces);
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

    std::atomic<int> num_threads(max_threads);
    split(root, aabbs, &num_threads);

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
