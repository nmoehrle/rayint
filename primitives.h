#include <limits>

struct AABB {
    math::Vec3f min;
    math::Vec3f max;
};

struct Tri {
    math::Vec3f a;
    math::Vec3f b;
    math::Vec3f c;
};

struct Ray {
    math::Vec3f origin;
    math::Vec3f dir;
    float tmin;
    float tmax;
};

inline
AABB operator+(AABB const & a, AABB const & b) {
    AABB aabb;
    for (std::size_t i = 0; i < 3; ++i) {
        aabb.min[i] = std::min(a.min[i], b.min[i]);
        aabb.max[i] = std::max(a.max[i], b.max[i]);
    }
    return aabb;
}

void calculate_aabb(Tri const & tri, AABB * aabb)  {
    for (std::size_t i = 0; i < 3; ++i) {
        aabb->min[i] = std::min(tri.a[i], std::min(tri.b[i], tri.c[i]));
        aabb->max[i] = std::max(tri.a[i], std::max(tri.b[i], tri.c[i]));
    }
}

float surface_area(AABB const & aabb) {
    float e0 = aabb.max[0] - aabb.min[0];
    float e1 = aabb.max[1] - aabb.min[1];
    float e2 = aabb.max[2] - aabb.min[2];
    return 2.0f * (e0 * e1 + e1 * e2 + e2 * e0);
}

float mid(AABB const & aabb, std::size_t d) {
    return (aabb.min[d] + aabb.max[d]) / 2.0f;
}

constexpr float inf = std::numeric_limits<float>::infinity();

bool intersect(Ray const & ray, AABB const & aabb) {
    float tmin = ray.tmin, tmax = ray.tmax;
    for (std::size_t i = 0; i < 3; ++i) {
        float t1 = (aabb.min[i] - ray.origin[i]) / ray.dir[i];
        float t2 = (aabb.max[i] - ray.origin[i]) / ray.dir[i];

        tmin = std::max(tmin, std::min(std::min(t1, t2), tmax));
        tmax = std::min(tmax, std::max(std::max(t1, t2), tmin));
    }
    return tmax > std::max(tmin, 0.0f);
}

bool intersect(Ray const & ray, Tri const & tri, float * t, math::Vec3f * bcoords = nullptr) {
    math::Vec3f v0 = tri.b - tri.a;
    math::Vec3f v1 = tri.c - tri.a;
    math::Vec3f normal = v1.cross(v0);
    if (normal.norm() < std::numeric_limits<float>::epsilon()) return false;
    normal.normalize();

    float cosine = normal.dot(ray.dir);
    if (std::abs(cosine) < std::numeric_limits<float>::epsilon()) return false;

    *t = -normal.dot(ray.origin - tri.a) / cosine;

    if (*t < ray.tmin || ray.tmax < *t) return false;
    math::Vec3f v2 = (ray.origin - tri.a) + *t * ray.dir;

    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;

    math::Vec3f uvw;
    uvw[0] = (d11 * d20 - d01 * d21) / denom;
    uvw[1] = (d00 * d21 - d01 * d20) / denom;
    uvw[2] = 1.0f - uvw[0] - uvw[1];

    if (bcoords != nullptr) {
        *bcoords = uvw;
    }

    return 0.0f <= uvw[0] && uvw[0] <= 1.0f
        && 0.0f <= uvw[1] && uvw[1] <= 1.0f
        && 0.0f <= uvw[2] && uvw[2] <= 1.0f;
}
