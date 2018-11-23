//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Gaspard Zoss, Alexandru Ichim
//
//   Copyright (C) 2016 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//   Edited 2017
//-----------------------------------------------------------------------------
#include "mesh_processing.h"
#include <Eigen/Dense>
#include <set>
#include <unordered_map>
#include <cmath>

namespace mesh_processing {

    using surface_mesh::Point;
    using surface_mesh::Scalar;
    using surface_mesh::Color;
    using surface_mesh::Vec2d;
    using surface_mesh::Vec3d;


    using std::min;
    using std::max;
    using std::cout;
    using std::endl;
    using std::pair;
    using std::vector;
    using std::unordered_map;


    MeshProcessing::MeshProcessing(const string &filename) {
        load_mesh(filename);
    }

// ======================================================================
// EXERCISE Begin
// ======================================================================

// ======================================================================
// EXERCISE 1.1 Mapping the Surface Boundary to 2D Circle
// ========================================================================
    void MeshProcessing::map_suface_boundary_to_circle() {

        Mesh::Vertex_property <Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0.0));
        int n_vertices = mesh_.n_vertices();
        //Homework starting from here

        // set everything to center of the circle
        for (const auto &v : mesh_.vertices()) v_texture[v] = Vec2d(.5, .5);

        /*
         * Idea: we construct a map that maps every boundary halfedge to the cumulated length along the boundary.
         * We start at some vertex on the boundary. Then we follow the halfedges on the boundary as long as we reach
         * the first vertex again. At the same time we sum the edge lengths and store the cumulated length of for each
         * boundary halfedge in the map.
         */
        Mesh::Vertex first_vertex;
        std::map<Mesh::Halfedge, double> boundary_edges;
        double cumulated_boundary_length = 0;

        // Find some vertex to start with, which is on the boundary...
        for (const auto &v : mesh_.vertices()) {
            if (mesh_.is_boundary(v)) {
                first_vertex = v;
                break;
            }
        }

        // ...then follow the boundary halfedges starting from it.
        // (if we traverse cw or ccw doesn't matter and depends on the mesh).
        Mesh::Halfedge h(-1);
        Mesh::Vertex v = first_vertex;
        do {
            bool found = false;

            // select next boundary halfedge not going backwards
            for (const auto &h1 : mesh_.halfedges(v)) {
                if (mesh_.is_boundary(h1) && mesh_.opposite_halfedge(h1) != h) {
                    h = h1;
                    found = true;
                    break;
                }
            }

            if (!found) throw std::runtime_error("Error traversing boundary vertices");

            // compute cumulated length and store the halfedge in the map
            cumulated_boundary_length += mesh_.edge_length(mesh_.edge(h));
            boundary_edges[h] = cumulated_boundary_length;
            v = mesh_.to_vertex(h);
        } while (v != first_vertex);

        const double total_length = cumulated_boundary_length;
        cout << "Found " << boundary_edges.size() << " boundary edges with a length of " << total_length << "." << endl;

        // Set vertex position on unit circle. We compute the angle at which the vertex appears on the unit circle and
        // map it to 2D coordinates using cosine and sine of that angle
        for (const auto& [halfedge, cumulated_length] : boundary_edges ) {
            double phi = 2*M_PI*cumulated_length/total_length;
            v_texture[mesh_.to_vertex(halfedge)] = (Vec2d(cos(phi), sin(phi)) / 2) + Vec2d(.5,.5);
        }
        //Homework stopping from here

        //Update the texture matrix￼
        texture_ = Eigen::MatrixXf(2, n_vertices);
        int j = 0;
        for (auto v: mesh_.vertices()) {
            texture_.col(j) << v_texture[v][0], v_texture[v][1];
            j++;
        }
    }

// ======================================================================
// EXERCISE 1.2 Iterative Solve Textures
// ========================================================================
    void MeshProcessing::iterative_solve_textures(int item_times) {
        int n_vertices = mesh_.n_vertices();
        Mesh::Vertex_property <Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0.5, 0.5));

        //Homework starting from here
        Mesh::Vertex_property <Vec2d> v_texture_new = mesh_.vertex_property<Vec2d>("v:texturenew", Vec2d(0, 0));

        Mesh::Edge_property <Scalar> e_weight = mesh_.edge_property<Scalar>("e:weight");

        for (const auto &vi : mesh_.vertices()) {
            float sum_edge_weights = 0.f;
            for (const auto &h : mesh_.halfedges(vi)) {
                const float h_weight = e_weight[mesh_.edge(h)];
                sum_edge_weights += h_weight;
                v_texture_new[vi] += h_weight * v_texture[mesh_.to_vertex(h)];
            }
            v_texture_new[vi] /= sum_edge_weights;
        }

        // overwrite texture with new texture
        for (const auto &v : mesh_.vertices()) {
            if (!mesh_.is_boundary(v)) {
                v_texture[v] = v_texture_new[v];
            }
        }

        cout << "Iterative Solve: Done" << endl;

        // clean up
        mesh_.remove_vertex_property(v_texture_new);


        //Homework stopping from here
        //Update the texture matrix
        texture_ = Eigen::MatrixXf(2, n_vertices);
        int j = 0;
        for (auto v: mesh_.vertices()) {
            texture_.col(j) << v_texture[v][0], v_texture[v][1];
            j++;
        }
    }

// ======================================================================
// EXERCISE 1.3 Direct Solve Textures
// ========================================================================
    void MeshProcessing::direct_solve_textures() {
        Mesh::Vertex_property <Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0.0));
        int n_vertices = mesh_.n_vertices();
        //Homework starting from here

        //Homework stopping from here
        //Update the texture matrix
        texture_ = Eigen::MatrixXf(2, n_vertices);
        int j = 0;
        for (auto v: mesh_.vertices()) {
            texture_.col(j) << v_texture[v][0], v_texture[v][1];
            j++;
        }
    }

// ======================================================================
// EXERCISE 2 Minimal Surfaces
// ======================================================================
    void MeshProcessing::minimal_surface() {

    }

// ======================================================================
// EXERCISE End
// ======================================================================


    void MeshProcessing::calc_weights() {
        calc_edges_weights();
        calc_vertices_weights();
    }

    void MeshProcessing::calc_uniform_mean_curvature() {
        Mesh::Vertex_property <Scalar> v_unicurvature =
                mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    }

    void MeshProcessing::calc_mean_curvature() {
        Mesh::Vertex_property <Scalar> v_curvature =
                mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
        Mesh::Edge_property <Scalar> e_weight =
                mesh_.edge_property<Scalar>("e:weight", 0.0f);
        Mesh::Vertex_property <Scalar> v_weight =
                mesh_.vertex_property<Scalar>("v:weight", 0.0f);
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    }

    void MeshProcessing::calc_gauss_curvature() {
        Mesh::Vertex_property <Scalar> v_gauss_curvature =
                mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
        Mesh::Vertex_property <Scalar> v_weight =
                mesh_.vertex_property<Scalar>("v:weight", 0.0f);
        // ------------- COPY YOUR FUNCTION FROM EXERCISE 4 ---------
    }

    void MeshProcessing::calc_edges_weights() {
        auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
        auto points = mesh_.vertex_property<Point>("v:point");

        Mesh::Halfedge h0, h1, h2;
        Point p0, p1, p2, d0, d1;

        for (auto e: mesh_.edges()) {
            double w = 0;
            e_weight[e] = 0.0;

            h0 = mesh_.halfedge(e, 0);
            p0 = points[mesh_.to_vertex(h0)];

            h1 = mesh_.halfedge(e, 1);
            p1 = points[mesh_.to_vertex(h1)];

            if (!mesh_.is_boundary(h0)) {
                h2 = mesh_.next_halfedge(h0);
                p2 = points[mesh_.to_vertex(h2)];
                d0 = p0 - p2;
                d1 = p1 - p2;
                w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));
            }

            if (!mesh_.is_boundary(h1)) {
                h2 = mesh_.next_halfedge(h1);
                p2 = points[mesh_.to_vertex(h2)];
                d0 = p0 - p2;
                d1 = p1 - p2;
                w += 1.0 / tan(acos(std::min(0.99f, std::max(-0.99f, dot(d0, d1)))));
            }

            w = w < 0 ? 0 : w;
            e_weight[e] = w * 0.5;
        }
    }

    void MeshProcessing::calc_vertices_weights() {
        Mesh::Face_around_vertex_circulator vf_c, vf_end;
        Mesh::Vertex_around_face_circulator fv_c;
        Scalar area;
        auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

        for (auto v: mesh_.vertices()) {
            area = 0.0;
            vf_c = mesh_.faces(v);

            if (!vf_c) {
                continue;
            }

            vf_end = vf_c;

            do {
                fv_c = mesh_.vertices(*vf_c);

                const Point &P = mesh_.position(*fv_c);
                ++fv_c;
                const Point &Q = mesh_.position(*fv_c);
                ++fv_c;
                const Point &R = mesh_.position(*fv_c);

                area += norm(cross(Q - P, R - P)) * 0.5f * 0.3333f;

            } while (++vf_c != vf_end);

            v_weight[v] = 0.5 / area;
        }
    }

    void MeshProcessing::load_mesh(const string &filename) {
        if (!mesh_.read(filename)) {
            std::cerr << "Mesh not found, exiting." << std::endl;
            exit(-1);
        }

        cout << "Mesh " << filename << " loaded." << endl;
        cout << "# of vertices : " << mesh_.n_vertices() << endl;
        cout << "# of faces : " << mesh_.n_faces() << endl;
        cout << "# of edges : " << mesh_.n_edges() << endl;

        // Compute the center of the mesh
        mesh_center_ = Point(0.0f, 0.0f, 0.0f);
        for (auto v: mesh_.vertices()) {
            mesh_center_ += mesh_.position(v);
        }
        mesh_center_ /= mesh_.n_vertices();

        // Compute the maximum distance from all points in the mesh and the center
        dist_max_ = 0.0f;
        for (auto v: mesh_.vertices()) {
            if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
                dist_max_ = distance(mesh_center_, mesh_.position(v));
            }
        }

        compute_mesh_properties();

        init_textures();

        // Store the original mesh, this might be useful for some computations
        mesh_init_ = mesh_;
    }

    void MeshProcessing::init_textures() {
        Mesh::Vertex_property <Vec2d> v_texture = mesh_.vertex_property<Vec2d>("v:texture", Vec2d(0, 0));
        int n_vertices = mesh_.n_vertices();
        texture_ = Eigen::MatrixXf(2, n_vertices);
        int j = 0;

        double min[3] = {1e10, 1e10, 1e10};
        double max[3] = {-1e10, -1e10, -1e10};
        for (auto v : mesh_.vertices()) {
            Point p = mesh_.position(v);
            for (int kd = 0; kd < 3; kd++) {
                if (p[kd] < min[kd]) min[kd] = p[kd];
                if (p[kd] > max[kd]) max[kd] = p[kd];
            }
        }
        for (auto v: mesh_.vertices()) {
            Point p = mesh_.position(v);
            v_texture[v][0] = (p[0] - min[0]) / (max[0] - min[0]);
            v_texture[v][1] = (p[1] - min[1]) / (max[1] - min[1]);
            texture_.col(j) << v_texture[v][0],
                    v_texture[v][1];
            j++;
        }
    }

    void MeshProcessing::compute_mesh_properties() {
        Mesh::Vertex_property <Point> vertex_normal =
                mesh_.vertex_property<Point>("v:normal");
        mesh_.update_face_normals();
        mesh_.update_vertex_normals();
        Mesh::Vertex_property <Color> v_color_valence =
                mesh_.vertex_property<Color>("v:color_valence",
                                             Color(1.0f, 1.0f, 1.0f));
        Mesh::Vertex_property <Color> v_color_unicurvature =
                mesh_.vertex_property<Color>("v:color_unicurvature",
                                             Color(1.0f, 1.0f, 1.0f));
        Mesh::Vertex_property <Color> v_color_curvature =
                mesh_.vertex_property<Color>("v:color_curvature",
                                             Color(1.0f, 1.0f, 1.0f));
        Mesh::Vertex_property <Color> v_color_gaussian_curv =
                mesh_.vertex_property<Color>("v:color_gaussian_curv",
                                             Color(1.0f, 1.0f, 1.0f));

        Mesh::Vertex_property <Scalar> vertex_valence =
                mesh_.vertex_property<Scalar>("v:valence", 0.0f);
        for (auto v: mesh_.vertices()) {
            vertex_valence[v] = mesh_.valence(v);
        }

        Mesh::Vertex_property <Scalar> v_unicurvature =
                mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
        Mesh::Vertex_property <Scalar> v_curvature =
                mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
        Mesh::Vertex_property <Scalar> v_gauss_curvature =
                mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);

        calc_weights();
        calc_uniform_mean_curvature();
        calc_mean_curvature();
        calc_gauss_curvature();
        color_coding(vertex_valence, &mesh_, v_color_valence, 100 /* bound */);
        color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
        color_coding(v_curvature, &mesh_, v_color_curvature);
        color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

        // get the mesh attributes and upload them to the GPU
        int j = 0;
        unsigned int n_vertices(mesh_.n_vertices());

        // Create big matrices to send the data to the GPU with the required
        // format
        color_valence_ = Eigen::MatrixXf(3, n_vertices);
        color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
        color_curvature_ = Eigen::MatrixXf(3, n_vertices);
        color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
        normals_ = Eigen::MatrixXf(3, n_vertices);
        points_ = Eigen::MatrixXf(3, n_vertices);
        indices_ = MatrixXu(3, mesh_.n_faces());

        for (auto f: mesh_.faces()) {
            std::vector<float> vv(3);
            int k = 0;
            for (auto v: mesh_.vertices(f)) {
                vv[k] = v.idx();
                ++k;
            }
            indices_.col(j) << vv[0], vv[1], vv[2];
            ++j;
        }

        j = 0;
        for (auto v: mesh_.vertices()) {
            points_.col(j) << mesh_.position(v).x,
                    mesh_.position(v).y,
                    mesh_.position(v).z;

            normals_.col(j) << vertex_normal[v].x,
                    vertex_normal[v].y,
                    vertex_normal[v].z;

            color_valence_.col(j) << v_color_valence[v].x,
                    v_color_valence[v].y,
                    v_color_valence[v].z;

            color_unicurvature_.col(j) << v_color_unicurvature[v].x,
                    v_color_unicurvature[v].y,
                    v_color_unicurvature[v].z;

            color_curvature_.col(j) << v_color_curvature[v].x,
                    v_color_curvature[v].y,
                    v_color_curvature[v].z;

            color_gaussian_curv_.col(j) << v_color_gaussian_curv[v].x,
                    v_color_gaussian_curv[v].y,
                    v_color_gaussian_curv[v].z;
            ++j;
        }
    }

    void MeshProcessing::color_coding(Mesh::Vertex_property <Scalar> prop, Mesh *mesh,
                                      Mesh::Vertex_property <Color> color_prop, int bound) {
        // Get the value array
        std::vector<Scalar> values = prop.vector();

        // discard upper and lower bound
        unsigned int n = values.size() - 1;
        unsigned int i = n / bound;
        std::sort(values.begin(), values.end());
        Scalar min_value = values[i], max_value = values[n - 1 - i];

        // map values to colors
        for (auto v: mesh->vertices()) {
            set_color(v, value_to_color(prop[v], min_value, max_value), color_prop);
        }
    }

    void MeshProcessing::set_color(Mesh::Vertex v, const Color &col,
                                   Mesh::Vertex_property <Color> color_prop) {
        color_prop[v] = col;
    }

    Color MeshProcessing::value_to_color(Scalar value, Scalar min_value, Scalar max_value) {
        Scalar v0, v1, v2, v3, v4;
        v0 = min_value + 0.0 / 4.0 * (max_value - min_value);
        v1 = min_value + 1.0 / 4.0 * (max_value - min_value);
        v2 = min_value + 2.0 / 4.0 * (max_value - min_value);
        v3 = min_value + 3.0 / 4.0 * (max_value - min_value);
        v4 = min_value + 4.0 / 4.0 * (max_value - min_value);

        Color col(1.0f, 1.0f, 1.0f);

        if (value < v0) {
            col = Color(0, 0, 1);
        } else if (value > v4) {
            col = Color(1, 0, 0);
        } else if (value <= v2) {
            if (value <= v1) { // [v0, v1]
                Scalar u = (value - v0) / (v1 - v0);
                col = Color(0, u, 1);
            } else { // ]v1, v2]
                Scalar u = (value - v1) / (v2 - v1);
                col = Color(0, 1, 1 - u);
            }
        } else {
            if (value <= v3) { // ]v2, v3]
                Scalar u = (value - v2) / (v3 - v2);
                col = Color(u, 1, 0);
            } else { // ]v3, v4]
                Scalar u = (value - v3) / (v4 - v3);
                col = Color(1, 1 - u, 0);
            }
        }
        return col;
    }

    MeshProcessing::~MeshProcessing() {}
}
