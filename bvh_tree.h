#include <queue>
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
    std::vector<std::size_t> const & faces;
    std::vector<math::Vec3f> const & vertices;
    std::vector<std::size_t> indices;

    Node *root;

public:
    ~BVHTree() {
        std::queue<Node*> q;
        q.push(root);
        while (!q.empty()) {
            Node *node = q.front(); q.pop();
            if (node == nullptr) continue;

            q.push(node->left);
            q.push(node->right);
            delete node;
        }
    };

    BVHTree(std::vector<std::size_t> const & faces, std::vector<math::Vec3f> const & vertices)
        : faces(faces), vertices(vertices) {
        std::vector<AABB> aabbs(faces.size() / 3);
        std::vector<AABB> right_aabbs(faces.size() / 3);
        for (std::size_t i = 0; i < aabbs.size(); ++i) {
            Tri tri;
            tri.a = vertices[faces[i * 3 + 0]];
            tri.b = vertices[faces[i * 3 + 1]];
            tri.c = vertices[faces[i * 3 + 2]];

            calculate_aabb(tri, &aabbs[i]);
        }
        indices.resize(aabbs.size());
        for (std::size_t i = 0; i < indices.size(); ++i) {
            indices[i] = i;
        }
        root = new Node(0, indices.size());

        std::queue<Node*> q;
        q.push(root);
        while (!q.empty()) {
            Node *node = q.front(); q.pop();
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

            if (min_cost < n && n > 5) {
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
                q.push(node->left);
                q.push(node->right);
            }
        }
    }

    struct Hit {
        float t;
        std::size_t idx;
        math::Vec3f bcoords;
    };

    void intersect(Ray const & ray, Node const * node, std::vector<Hit> * hits) {
        for (std::size_t i = node->first; i < node->last; ++i) {
            Tri tri;
            std::size_t face_idx = indices[i] * 3;
            tri.a = vertices[faces[face_idx + 0]];
            tri.b = vertices[faces[face_idx + 1]];
            tri.c = vertices[faces[face_idx + 2]];

            Hit hit;
            if (::intersect(ray, tri, &hit.t, &hit.bcoords)) {
                hit.idx = indices[i];
                hits->push_back(hit);
            }
        }
    }

    bool intersect(Ray const & ray, Hit * hit) {
        std::vector<Hit> hits;
        std::stack<Node const *> s;

        s.push(root);
        while (!s.empty()) {
            Node const *node = s.top(); s.pop();
            if (node == nullptr) continue;
            if (::intersect(ray, node->aabb)) {
                if (node->is_leaf()) {
                    intersect(ray, node, &hits);
                } else {
                    s.push(node->left);
                    s.push(node->right);
                }
            }
        }

        if (hits.empty()) return false;

        std::sort(hits.begin(), hits.end(),
            [] (Hit const & first, Hit const & second) -> bool {
                return first.t < second.t;
            }
        );
        *hit = hits[0];
        return true;
    }
};
