#include <deque>
#include <stack>
#include <algorithm>

#include <math/vector.h>

#include "primitives.h"

class BVHTree {
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

public:
    ~BVHTree() {
        std::deque<Node*> q;
        q.push_back(root);
        while (!q.empty()) {
            Node *node = q.back(); q.pop_back();
            if (node == nullptr) continue;

            q.push_back(node->left);
            q.push_back(node->right);
            delete node;
        }
    };

    BVHTree(std::vector<std::size_t> const & faces, std::vector<math::Vec3f> const & vertices) {
        std::vector<AABB> aabbs(faces.size() / 3);
        std::vector<AABB> right_aabbs(faces.size() / 3);
        std::vector<Tri> ttris(faces.size() / 3);
        for (std::size_t i = 0; i < aabbs.size(); ++i) {
            ttris[i].a = vertices[faces[i * 3 + 0]];
            ttris[i].b = vertices[faces[i * 3 + 1]];
            ttris[i].c = vertices[faces[i * 3 + 2]];

            calculate_aabb(ttris[i], &aabbs[i]);
        }
        indices.resize(aabbs.size());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }
        root = new Node(0, indices.size());

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
                            || (mid(aabbs[first], d) == mid(aabbs[second], d) && first < second);
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
                    if (cost < min_cost) {
                        min_cost = cost;
                        split = std::make_pair(d, i);
                    }

                    left_aabb = left_aabb + aabbs[indices[i]];
                }
            }

            if (min_cost < n) {
                std::size_t d, i;
                std::tie(d, i) = split;
                std::sort(&indices[node->first], &indices[node->last],
                    [&aabbs, d] (std::size_t first, std::size_t second) -> bool {
                        return mid(aabbs[first], d) < mid(aabbs[second], d)
                            || (mid(aabbs[first], d) == mid(aabbs[second], d) && first < second);
                    }
                );

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

    struct Hit {
        float t;
        std::size_t idx;
        math::Vec3f bcoords;
    };

    bool intersect(Ray const & ray, Node const * node, Hit * hit) {
        bool ret = false;
        for (std::size_t i = node->first; i < node->last; ++i) {
            float t;
            math::Vec3f bcoords;
            if (::intersect(ray, tris[i], &t, &bcoords)) {
                if (t > hit->t) continue;
                hit->idx = indices[i];
                hit->t = t;
                hit->bcoords = bcoords;
                ret = true;
            }
        }
        return ret;
    }

    bool intersect(Ray ray, Hit * hit_ptr) {
        Hit hit;
        hit.t = std::numeric_limits<float>::infinity();
        std::stack<Node const *> s;

        s.push(root);
        while (!s.empty()) {
            Node const *node = s.top(); s.pop();
            if (node == nullptr) continue;
            if (::intersect(ray, node->aabb)) {
                if (node->is_leaf()) {
                    if (intersect(ray, node, &hit)) {
                        ray.tmax = hit.t;
                    }
                } else {
                    s.push(node->left);
                    s.push(node->right);
                }
            }
        }

        if (hit_ptr != nullptr) {
            *hit_ptr = hit;
        }

        return hit.t < std::numeric_limits<float>::infinity();
    }
};
