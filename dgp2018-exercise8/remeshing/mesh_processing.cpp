//=============================================================================
//
//   Code framework for the lecture
//
//   "Digital 3D Geometry Processing"
//
//   Copyright (C) 2017 by Computer Graphics and Geometry Laboratory,
//         EPF Lausanne
//
//-----------------------------------------------------------------------------
#define _USE_MATH_DEFINES

#include "mesh_processing.h"
#include <set>
#include <cmath>

namespace mesh_processing {

    using surface_mesh::Point;
    using surface_mesh::Scalar;
    using surface_mesh::Color;
    using std::min;
    using std::max;
    using std::cout;
    using std::endl;

    MeshProcessing::MeshProcessing(const string &filename) {
        load_mesh(filename);
    }

    MeshProcessing::~MeshProcessing() {
        // TODO
    }

    void MeshProcessing::remesh(const REMESHING_TYPE &remeshing_type,
                                const int &num_iterations) {
        calc_weights();
        calc_mean_curvature();
        calc_uniform_mean_curvature();
        calc_gauss_curvature();
        calc_target_length(remeshing_type);

        // main remeshing loop
        for (int i = 0; i < num_iterations; ++i) {
            split_long_edges();
//            collapse_short_edges();
//            equalize_valences();
//            tangential_relaxation();
        }
    }

    void MeshProcessing::calc_target_length(const REMESHING_TYPE &remeshing_type) {
        Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
        Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
        Scalar length;
        Scalar mean_length;
        Scalar H;
        Scalar K;

        Mesh::Vertex_property <Scalar> curvature = mesh_.vertex_property<Scalar>("v:meancurvature", 0);
        Mesh::Vertex_property <Scalar> gauss_curvature = mesh_.vertex_property<Scalar>("v:gausscurvature", 0);
        Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
        Mesh::Vertex_property <Scalar> target_new_length = mesh_.vertex_property<Scalar>("v:newlength", 0);

        // our helpers
        Mesh::Vertex_property <Scalar> max_curvature = mesh_.vertex_property<Scalar>("v:maxcurvature", 0);
        Mesh::Vertex_property <Scalar> v_new_target_length = mesh_.vertex_property<Scalar>("v:vnewtargetlength", 0);

        // user specified target length
        const float TARGET_LENGTH = 2.0;

        if (remeshing_type == AVERAGE) {
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
                target_length[*v_it] = TARGET_LENGTH;

        } else if (remeshing_type == CURV) {
            // ------------- IMPLEMENT HERE ---------
            // Get the maximal curvature at each vertex (use the precomputed mean and gaussian curvature)
            // Calculate the desired edge length as the TARGET_LENGTH divided by the maximal curvature at each vertex, and assign it to the property target_length
            // Smooth the maximal curvature uniformly, use the property vnewtargetlength_ to store the smoothed values intermediately
            // Rescale the property target_new_length such that it's mean equals the user specified TARGET_LENGTH
            // ------------- IMPLEMENT HERE ---------

            // calculate desired length
            for (const auto &v : mesh_.vertices()) {
                length = 1.0;
                if (!mesh_.is_boundary(v)) {
                    max_curvature[v] =
                            curvature[v] + sqrtf(curvature[v] * curvature[v] - gauss_curvature[v]);
                    target_length[v] = TARGET_LENGTH / max_curvature[v];
                }
                target_length[v] = length;
            }

            // smooth desired length uniformly
            for (int i = 0; i < 5; i++) {
                for (const auto &v1 : mesh_.vertices()) {
                    v_new_target_length[v1] = 0.;

                    int n = 0; // number of neighbors
                    for (const auto &v2 : mesh_.vertices(v1)) {
                        v_new_target_length[v1] += target_length[v2] - target_length[v1];
                        n += 1;
                    }
                    v_new_target_length[v1] /= n;
                }
            }

            // calculate mean of v_new_target_length
            Scalar mean_length(0.);
            unsigned int N = mesh_.n_vertices();
            for (const auto &v : mesh_.vertices()) {
                mean_length += v_new_target_length[v] / N;
            }

            // rescale desired length
            for (const auto &v : mesh_.vertices()) {
                target_new_length[v] = v_new_target_length[v] * (TARGET_LENGTH / mean_length);
            }
        }
    }

    void MeshProcessing::split_long_edges() {
        Mesh::Vertex v0, v1, v;
        bool finished;
        int i;

        const float upper_ratio = 4 / 3.f;

        Mesh::Vertex_property <Point> normals = mesh_.vertex_property<Point>("v:normal");
        Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

        for (finished = false, i = 0; !finished && i < 100; ++i) {
            finished = true;
            int c = 0;
            // ------------- IMPLEMENT HERE ---------
            // INSERT CODE:
            //  Compute the desired length as the mean between the property target_length of two vertices of the edge
            //  If the edge is longer than 4/3 * desired length
            //		add the midpoint to the mesh
            //		set the interpolated normal and interpolated vtargetlength_ property to the vertex
            //		split the edge with this vertex (use openMesh function split)
            // Leave the loop running until no splits are done (use the finished variable)
            // ------------- IMPLEMENT HERE ---------

            for (const auto &edge : mesh_.edges()) {
                v0 = mesh_.vertex(edge, 0);
                v1 = mesh_.vertex(edge, 1);
                float desired_length = .5f * target_length[v0] + .5f * target_length[v1];
                if (mesh_.edge_length(edge) > upper_ratio * desired_length) {
                    finished = false;
                    v = mesh_.split(edge,
                                    mesh_.position(v0) + .5 * (mesh_.position(v1) - mesh_.position(v0)));
                    normals[v] = .5f * normals[v0] + .5f * normals[v1];
                    target_length[v] = desired_length;
                    c++;
                }
            }
            cout << "Splitted " << c << " long edges in run " << i << "." << endl;
            mesh_.update_vertex_normals();
        }
    }

    void MeshProcessing::collapse_short_edges() {
        Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
        Mesh::Vertex v0, v1;
        Mesh::Halfedge h01, h10;
        bool finished, b0, b1;
        int i;
        bool hcol01, hcol10;

        Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

        for (finished = false, i = 0; !finished && i < 100; ++i) {
            finished = true;

            for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
                if (!mesh_.is_deleted(*e_it)) // might already be deleted
                {
                    // ------------- IMPLEMENT HERE ---------
                    // INSERT CODE:
                    // Compute the desired length as the mean between the property vtargetlength_ of two vertices of the edge
                    // If the edge is shorter than 4/5 of the desired length
                    //		Check if halfedge connects a boundary vertex with a non-boundary vertex. If so, don't collapse.
                    //		Check if halfedges collapsible
                    //		Select the halfedge to be collapsed if at least one halfedge can be collapsed
                    //		Collapse the halfedge
                    // Leave the loop running until no collapse has been done (use the finished variable)
                    // ------------- IMPLEMENT HERE ---------
                }
            }
        }

        mesh_.garbage_collection();

        if (i == 100) std::cerr << "collapse break\n";
    }

    void MeshProcessing::equalize_valences() {
        Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
        Mesh::Vertex v0, v1, v2, v3;
        Mesh::Halfedge h;
        int val0, val1, val2, val3;
        int val_opt0, val_opt1, val_opt2, val_opt3;
        int ve0, ve1, ve2, ve3, ve_before, ve_after;
        bool finished;
        int i;


        // flip all edges
        for (finished = false, i = 0; !finished && i < 100; ++i) {
            finished = true;

            for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
                if (!mesh_.is_boundary(*e_it)) {
                    // ------------- IMPLEMENT HERE ---------
                    //  Extract valences of the four vertices involved to an eventual flip.
                    //  Compute the sum of the squared valence deviances before flip
                    //  Compute the sum of the squared valence deviances after an eventual flip
                    //  If valence deviance is decreased and flip is possible, flip the vertex
                    //  Leave the loop running until no collapse has been done (use the finished variable)
                    // ------------- IMPLEMENT HERE ---------
                }
            }
        }

        if (i == 100) std::cerr << "flip break\n";
    }

    void MeshProcessing::tangential_relaxation() {
        Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
        Mesh::Vertex_around_vertex_circulator vv_c, vv_end;
        int valence;
        Point u, n;
        Point laplace;

        Mesh::Vertex_property <Point> normals = mesh_.vertex_property<Point>("v:normal");
        Mesh::Vertex_property <Point> update = mesh_.vertex_property<Point>("v:update");

        // smooth
        for (int iters = 0; iters < 10; ++iters) {
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
                if (!mesh_.is_boundary(*v_it)) {
                    // ------------- IMPLEMENT HERE ---------
                    //  Compute uniform laplacian curvature approximation vector
                    //  Compute the tangential component of the laplacian vector and move the vertex
                    //  Store smoothed vertex location in the update vertex property.
                    // ------------- IMPLEMENT HERE ---------
                }
            }

            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
                if (!mesh_.is_boundary(*v_it))
                    mesh_.position(*v_it) += update[*v_it];
        }
    }

    void MeshProcessing::calc_uniform_mean_curvature() {
        Mesh::Vertex_property <Scalar> v_unicurvature =
                mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, approximate mean curvature using
        // the length of the uniform Laplacian approximation
        // Save your approximation in unicurvature vertex property of the mesh.
        // ------------- IMPLEMENT HERE ---------
    }

    void MeshProcessing::calc_mean_curvature() {
        Mesh::Vertex_property <Scalar> v_curvature =
                mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
        Mesh::Edge_property <Scalar> e_weight =
                mesh_.edge_property<Scalar>("e:weight", 0.0f);
        Mesh::Vertex_property <Scalar> v_weight =
                mesh_.vertex_property<Scalar>("v:weight", 0.0f);

        // ------------- IMPLEMENT HERE ---------
        // For all non-boundary vertices, approximate the mean curvature using
        // the length of the Laplace-Beltrami approximation.
        // Save your approximation in v_curvature vertex property of the mesh.
        // Use the weights from calc_weights(): e_weight and v_weight
        // ------------- IMPLEMENT HERE ---------
    }

    void MeshProcessing::calc_gauss_curvature() {
        Mesh::Vertex_property <Scalar> v_gauss_curvature =
                mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
        Mesh::Vertex_property <Scalar> v_weight =
                mesh_.vertex_property<Scalar>("v:weight", 0.0f);

        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, approximate Gaussian curvature,
        // and store it in the vertex property v_gauss_curvature.
        // Hint: When calculating angles out of cross products make sure the value
        // you pass to the acos function is between -1.0 and 1.0.
        // Use the v_weight property for the area weight.
        // ------------- IMPLEMENT HERE ---------
    }

    void MeshProcessing::calc_weights() {
        calc_edges_weights();
        calc_vertices_weights();
    }

    void MeshProcessing::calc_edges_weights() {
        auto e_weight = mesh_.edge_property<Scalar>("e:weight", 0.0f);
        auto points = mesh_.vertex_property<Point>("v:point");

        Mesh::Halfedge h0, h1, h2;
        Point p0, p1, p2, d0, d1;

        for (auto e : mesh_.edges()) {
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
                e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
            }

            if (!mesh_.is_boundary(h1)) {
                h2 = mesh_.next_halfedge(h1);
                p2 = points[mesh_.to_vertex(h2)];
                d0 = p0 - p2;
                d1 = p1 - p2;
                e_weight[e] += dot(d0, d1) / norm(cross(d0, d1));
            }
        }
    }

    void MeshProcessing::calc_vertices_weights() {
        Mesh::Face_around_vertex_circulator vf_c, vf_end;
        Mesh::Vertex_around_face_circulator fv_c;
        Scalar area;
        auto v_weight = mesh_.vertex_property<Scalar>("v:weight", 0.0f);

        for (auto v : mesh_.vertices()) {
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
        for (auto v : mesh_.vertices()) {
            mesh_center_ += mesh_.position(v);
        }
        mesh_center_ /= mesh_.n_vertices();

        // Compute the maximum distance from all points in the mesh and the center
        dist_max_ = 0.0f;
        for (auto v : mesh_.vertices()) {
            if (distance(mesh_center_, mesh_.position(v)) > dist_max_) {
                dist_max_ = distance(mesh_center_, mesh_.position(v));
            }
        }

        compute_mesh_properties();

        // Store the original mesh, this might be useful for some computations
        mesh_init_ = mesh_;
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
        for (auto v : mesh_.vertices()) {
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
        color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
                     8 /* max */);
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

        for (auto f : mesh_.faces()) {
            std::vector<float> vv(3);
            int k = 0;
            for (auto v : mesh_.vertices(f)) {
                vv[k] = v.idx();
                ++k;
            }
            indices_.col(j) << vv[0], vv[1], vv[2];
            ++j;
        }

        j = 0;
        for (auto v : mesh_.vertices()) {
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
                                      Mesh::Vertex_property <Color> color_prop, Scalar min_value,
                                      Scalar max_value, int bound) {
        // Get the value array
        std::vector<Scalar> values = prop.vector();

        if (min_value == 0.0 && max_value == 0.0) {
            // discard upper and lower bound
            unsigned int n = values.size() - 1;
            unsigned int i = n / bound;
            std::sort(values.begin(), values.end());
            min_value = values[i];
            max_value = values[n - 1 - i];
        }

        // map values to colors
        for (auto v : mesh->vertices()) {
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


}


