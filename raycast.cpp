#include <iostream>

#include <util/timer.h>
#include <util/arguments.h>
#include <mve/mesh_io_ply.h>
#include <mve/image_io.h>

#include "bvh_tree.h"

struct Arguments {
    std::string in_view;
    std::string image;
    std::string in_mesh;
    std::string out_image;
};

Arguments parse_args(int argc, char **argv) {
    util::Arguments args;
    args.set_exit_on_error(true);
    args.set_nonopt_maxnum(4);
    args.set_nonopt_minnum(4);
    args.set_usage("Usage: " + std::string(argv[0]) + " [OPTS] IN_VIEW IMAGE IN_MESH OUT_IMAGE");
    args.set_description("TODO");
    args.parse(argc, argv);

    Arguments conf;
    conf.in_view = args.get_nth_nonopt(0);
    conf.image = args.get_nth_nonopt(1);
    conf.in_mesh = args.get_nth_nonopt(2);
    conf.out_image = args.get_nth_nonopt(3);

    for (util::ArgResult const* i = args.next_option();
         i != 0; i = args.next_option()) {
        switch (i->opt->sopt) {

        default:
            throw std::invalid_argument("Invalid option");
        }
    }
    return conf;
}
int main(int argc, char **argv) {
    Arguments conf = parse_args(argc, argv);

    mve::TriangleMesh::Ptr mesh;
    try {
        mesh = mve::geom::load_ply_mesh(conf.in_mesh);
    } catch (std::exception& e) {
        std::cerr << "\tCould not load mesh: "<< e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    mve::View::Ptr view;
    try {
        view = mve::View::create(conf.in_view);
    } catch (std::exception& e) {
        std::cerr << "\tCould not load view: "<< e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (!view->has_image(conf.image)) {
        std::cerr << "\tView has no image: " << conf.image << std::endl;
        std::exit(EXIT_FAILURE);
    }
    mve::View::ImageProxy const * proxy = view->get_image_proxy(conf.image);
    mve::CameraInfo const & camera = view->get_camera();

    std::vector<unsigned int> const & mfaces = mesh->get_faces();
    std::vector<std::size_t> faces(mfaces.begin(), mfaces.end());
    std::vector<math::Vec3f> const & vertices = mesh->get_vertices();
    std::vector<math::Vec4f> const & colors = mesh->get_vertex_colors();

    BVHTree bvhtree(faces, vertices);

    math::Vec3f origin;
    camera.fill_camera_pos(*origin);
    math::Matrix3f invproj;
    camera.fill_inverse_calibration(*invproj, proxy->width, proxy->height);
    math::Matrix3f c2w_rot;
    camera.fill_cam_to_world_rot(*c2w_rot);

    util::WallTimer timer;
    mve::ByteImage::Ptr image = mve::ByteImage::create(proxy->width, proxy->height, 3);
    #pragma omp parallel for
    for (int y = 0; y < image->width(); ++y) {
        for (int x = 0; x < image->width(); ++x) {
            Ray ray;
            ray.origin = origin;
            math::Vec3f v = invproj * math::Vec3f ((float)x + 0.5f, (float)y + 0.5f, 1.0f);
            ray.dir = c2w_rot.mult(v.normalized());
            ray.tmin = 0.0f;
            ray.tmax = std::numeric_limits<float>::infinity();

            BVHTree::Hit hit;
            if (bvhtree.intersect(ray, &hit)) {
                math::Vec4f const & c1 = colors[faces[hit.idx * 3 + 0]];
                math::Vec4f const & c2 = colors[faces[hit.idx * 3 + 1]];
                math::Vec4f const & c3 = colors[faces[hit.idx * 3 + 2]];
                math::Vec3f const & w = hit.bcoords;
                for (std::size_t c = 0; c < 3; ++c) {
                    image->at(x, y, c) = 255.0f *
                        math::interpolate(c1[c], c2[c], c3[c], w[0], w[1], w[2]);
                }
            }
        }
    }
    std::cout << "Raycasting took: " << timer.get_elapsed() << " ms" << std::endl;

    mve::image::save_png_file(image, conf.out_image);
    return EXIT_SUCCESS;
}
