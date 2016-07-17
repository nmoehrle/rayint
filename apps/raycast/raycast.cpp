/*
 * Copyright (C) 2015, Nils Moehrle
 * All rights reserved.
 *
 * This software may be modified and distributed under the terms
 * of the BSD 3-Clause license. See the LICENSE.txt file for details.
 */

#include <iostream>

#include <util/timer.h>
#include <util/arguments.h>
#include <mve/mesh_io_ply.h>
#include <mve/scene.h>
#include <future>

#include <acc/bvh_tree.h>

typedef unsigned int uint;

typedef acc::BVHTree<uint, math::Vec3f> BVHTree;

struct Arguments {
    std::string in_scene;
    std::string in_mesh;
    std::string out_prefix;

    uint width = 1920;
    uint height = 1080;

    bool normals;
    bool depth;
};

Arguments parse_args(int argc, char **argv) {
    util::Arguments args;
    args.set_exit_on_error(true);
    args.set_nonopt_maxnum(3);
    args.set_nonopt_minnum(3);
    args.set_usage("Usage: " + std::string(argv[0]) + " [OPTS] IN_SCENE IN_MESH OUT_PREFIX");
    args.set_description("TODO");
    //args.add_option('n', "normal-map", false, "write out normal map (OUT_PREFIX-normals) [false]");
    //args.add_option('d', "depth-map", false, "write out depth map (OUT_PREFIX-depth) [false]");
    args.parse(argc, argv);

    Arguments conf;
    conf.in_scene = args.get_nth_nonopt(0);
    conf.in_mesh = args.get_nth_nonopt(1);
    conf.out_prefix = args.get_nth_nonopt(2);
    conf.normals = false;
    conf.depth = false;

    for (util::ArgResult const* i = args.next_option();
         i != 0; i = args.next_option()) {
        switch (i->opt->sopt) {
        case 'n': conf.normals = true; break;
        case 'd': conf.depth = true; break;
        default:
            throw std::invalid_argument("Invalid option");
        }
    }
    return conf;
}
int main(int argc, char **argv) {
    Arguments args = parse_args(argc, argv);

    mve::TriangleMesh::Ptr mesh;
    try {
        mesh = mve::geom::load_ply_mesh(args.in_mesh);
    } catch (std::exception& e) {
        std::cerr << "\tCould not load mesh: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (!mesh->has_vertex_colors()) {
        std::cerr << "\tMesh has no vertex colors" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    mve::Scene::Ptr scene;
    try {
        scene = mve::Scene::create(args.in_scene);
    } catch (std::exception& e) {
        std::cerr << "Could not open scene: " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::vector<mve::View::Ptr> views;
    {
        mve::Scene::ViewList const & aviews = scene->get_views();
        views.reserve(aviews.size());
        for (mve::View::Ptr view : aviews) {
            if (view == nullptr) continue;
            if (!view->is_camera_valid()) continue;
            views.push_back(view);
        }
    }

    std::size_t num_views = views.size();

    if (num_views == 0) {
        std::cerr << "Scene does not contain valid views." << std::endl;
        std::exit(EXIT_FAILURE);
    }
    mesh->ensure_normals(false, true);
    std::vector<uint> const & faces = mesh->get_faces();
    std::vector<math::Vec3f> const & vertices = mesh->get_vertices();
    std::vector<math::Vec4f> const & colors = mesh->get_vertex_colors();
    std::vector<math::Vec3f> const & normals = mesh->get_vertex_normals();

    util::WallTimer timer;
    std::cout << "Building BVH from " << faces.size() / 3 << " faces... " << std::flush;
    BVHTree bvhtree(faces, vertices);
    std::cout << "done. (Took: " << timer.get_elapsed() << " ms)" << std::endl;

    std::vector<std::future<void> > futures;
    futures.reserve(num_views);

    timer.reset();
    std::cout << "Raycasting... " << std::flush;
    for (mve::View::Ptr view : views) {
        mve::CameraInfo const & camera = view->get_camera();

        math::Vec3f origin;
        camera.fill_camera_pos(*origin);
        math::Matrix3f invproj;
        camera.fill_inverse_calibration(*invproj, args.width, args.height);
        math::Matrix3f c2w_rot;
        camera.fill_cam_to_world_rot(*c2w_rot);

        mve::ByteImage::Ptr uvmap = mve::ByteImage::create(args.width, args.height, 3);
        mve::ByteImage::Ptr normalmap = mve::ByteImage::create(args.width, args.height, 3);
        mve::FloatImage::Ptr depthmap = mve::FloatImage::create(args.width, args.height, 1);

        #pragma omp parallel for
        for (int y = 0; y < uvmap->height(); ++y) {
            for (int x = 0; x < uvmap->width(); ++x) {
                BVHTree::Ray ray;
                ray.origin = origin;
                math::Vec3f v = invproj * math::Vec3f ((float)x + 0.5f, (float)y + 0.5f, 1.0f);
                ray.dir = c2w_rot.mult(v.normalized()).normalize();
                ray.tmin = 0.0f;
                ray.tmax = std::numeric_limits<float>::infinity();

                BVHTree::Hit hit;
                if (bvhtree.intersect(ray, &hit)) {
                    math::Vec3f const & n1 = normals[faces[hit.idx * 3 + 0]];
                    math::Vec3f const & n2 = normals[faces[hit.idx * 3 + 1]];
                    math::Vec3f const & n3 = normals[faces[hit.idx * 3 + 2]];
                    math::Vec4f const & c1 = colors[faces[hit.idx * 3 + 0]];
                    math::Vec4f const & c2 = colors[faces[hit.idx * 3 + 1]];
                    math::Vec4f const & c3 = colors[faces[hit.idx * 3 + 2]];
                    math::Vec3f const & w = hit.bcoords;
                    math::Vec4f color = math::interpolate(c1, c2, c3, w[0], w[1], w[2]);
                    math::Vec3f normal = math::interpolate(n1, n2, n3, w[0], w[1], w[2]).normalize();
                    for (std::size_t c = 0; c < 3; ++c) {
                        uvmap->at(x, y, c) = 255.0f * color[c];
                        normalmap->at(x, y, c) = 255.0f * (0.5f + normal[c] / 2.0f);
                    }
                    depthmap->at(x, y, 0) = (hit.t * ray.dir).norm();
                }
            }
        }

        view->set_image(uvmap, args.out_prefix + "-uv");
        view->set_image(normalmap, args.out_prefix + "-normals");
        view->set_image(depthmap, args.out_prefix + "-depth");
        futures.push_back(std::async(std::launch::async, [view] {
            view->save_view();
            view->cache_cleanup();
        }));
    }
    std::cout << "done. (Took: " << timer.get_elapsed() << " ms)" << std::endl;

    return EXIT_SUCCESS;
}
