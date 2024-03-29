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
#include <Eigen/Geometry>

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

    MeshProcessing::~MeshProcessing() = default;

    void MeshProcessing::remesh(const REMESHING_TYPE &remeshing_type,
                                const int &num_iterations) {
        calc_weights();
        calc_mean_curvature();
        calc_uniform_mean_curvature();
        calc_gauss_curvature();
        calc_max_curvature();
        calc_target_length(remeshing_type);

        // main remeshing loop
        for (int i = 0; i < num_iterations; ++i) {
            cout << "----------- iteration " << i + 1 << " out of " << num_iterations << " -----------" << endl;
            split_long_edges();
            collapse_short_edges();
            equalize_valences();
            tangential_relaxation();
        }
    }

    void MeshProcessing::calc_target_length(const REMESHING_TYPE &remeshing_type) {
        Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
        Scalar length;

        Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);
        Mesh::Vertex_property <Scalar> new_target_length = mesh_.vertex_property<Scalar>("v:new_length", 0);

        // we compute the max curvature at the same time as mean and gauss curvature are computed.
        Mesh::Vertex_property <Scalar> max_curvature = mesh_.vertex_property<Scalar>("v:max_curvature", 0);

        // user specified target length
        const float TARGET_LENGTH = 0.75;

        if (remeshing_type == AVERAGE) {
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
                target_length[*v_it] = TARGET_LENGTH;

        } else if (remeshing_type == CURV) {
            // ------------- IMPLEMENT HERE ---------
            // Get the maximal curvature at each vertex (use the precomputed mean and gaussian curvature)
            // Calculate the desired edge length as the TARGET_LENGTH divided by the maximal curvature at each vertex, and assign it to the property target_length
            // Smooth the maximal curvature uniformly, use the property vnewtargetlength_ to store the smoothed values intermediately
            // Rescale the property new_target_length such that it's mean equals the user specified TARGET_LENGTH
            // ------------- IMPLEMENT HERE ---------

            // calculate desired length
            for (const auto &v : mesh_.vertices()) {
                length = TARGET_LENGTH;
                if (mesh_.is_boundary(v)) {
                    length = TARGET_LENGTH;
                } else {
                    length = TARGET_LENGTH / max_curvature[v];
                }
                target_length[v] = length;
            }

            //uniform smoothing of desired length
            const int smoothing_iterations = 5;
            for (int i = 0; i < smoothing_iterations; i++) {
                // apply uniform smoothing on target_length
                for (const auto &vi : mesh_.vertices()) {
                    float local_target_length_avg = 0.f;
                    for (const auto &vj : mesh_.vertices(vi)) {
                        local_target_length_avg += target_length[vj] - target_length[vi];
                    }
                    local_target_length_avg /= mesh_.valence(vi);
                    new_target_length[vi] = target_length[vi] + local_target_length_avg;
                    assert(target_length[vi] > 0);
                }

                // copy the new target length to target length
                for (const auto &v : mesh_.vertices()) {
                    target_length[v] = new_target_length[v];
                }

            }

            // calculate mean of target_length
            float mean_length = 0;
            for (const auto &v : mesh_.vertices()) {
                mean_length += target_length[v];
            }
            mean_length /= mesh_.n_vertices();

            // rescale desired length
            const double rescale_factor = TARGET_LENGTH / mean_length;
            for (const auto &v : mesh_.vertices())
                target_length[v] *= rescale_factor;
        }

        // Compute color coding of target length (for debugging and visualization)
        Mesh::Vertex_property <Color> c = mesh_.get_vertex_property<Color>("v:color_target_length");
        color_coding(target_length, &mesh_, c);
        // and write it to the GPU matrix
        int j = 0;
        for (const auto &v : mesh_.vertices()) {
            color_target_length_.col(j++) << c[v].x, c[v].y, c[v].z;
        }
    }

    void MeshProcessing::split_long_edges() {
        Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
        Mesh::Vertex v0, v1, v;
        bool finished;
        int i;
        int c = 0;

        const double upper_ratio = 4. / 3;

        Mesh::Vertex_property <Point> normals = mesh_.vertex_property<Point>("v:normal");
        Mesh::Vertex_property <Scalar> target_length = mesh_.vertex_property<Scalar>("v:length", 0);

        for (finished = false, i = 0; !finished && i < 100; ++i) {
            finished = true;
            // ------------- IMPLEMENT HERE ---------
            // INSERT CODE:
            //  Compute the desired length as the mean between the property target_length of two vertices of the edge
            //  If the edge is longer than 4/3 * desired length
            //		add the midpoint to the mesh
            //		set the interpolated normal and interpolated vtargetlength_ property to the vertex
            //		split the edge with this vertex (use openMesh function split)
            // Leave the loop running until no splits are done (use the finished variable)
            // ------------- IMPLEMENT HERE ---------

            for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
                v0 = mesh_.vertex(*e_it, 0);
                v1 = mesh_.vertex(*e_it, 1);
                float desired_length = .5f * target_length[v0] + .5f * target_length[v1];
                if (mesh_.edge_length(*e_it) > upper_ratio * desired_length) {
                    finished = false;
                    v = mesh_.split(*e_it, mesh_.position(v0) + .5 * (mesh_.position(v1) - mesh_.position(v0)));
                    normals[v] = .5 * normals[v0] + .5 * normals[v1];
                    target_length[v] = desired_length;
                    c++; // counts the number of splitted edges
                }
            }
        }
        cout << "Split " << c << " long edges in " << i << " iterations." << endl;
    }

    void MeshProcessing::collapse_short_edges() {
        Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
        Mesh::Vertex v0, v1;
        Mesh::Halfedge h01, h10, h;
        bool finished;
        int i, c = 0;

        const double lower_ratio = 4. / 5;

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
                    h01 = mesh_.halfedge(*e_it, 0);
                    h10 = mesh_.halfedge(*e_it, 1);

                    // to vertex give the vertex the halfedge points to
                    v0 = mesh_.to_vertex(h10);
                    v1 = mesh_.to_vertex(h01);

                    float desired_length = .5f * target_length[v0] + .5f * target_length[v1];

                    if (mesh_.edge_length(*e_it) < lower_ratio * desired_length) {

                        h = Mesh::Halfedge(); // sanity check, make h invalid before the selection

                        if (mesh_.is_boundary(v0) && !mesh_.is_boundary(v1)) { // collapse towards the boundary
                            h = h10;
                        } else if (!mesh_.is_boundary(v0) && mesh_.is_boundary(v1)) { // collapse towards the boundary
                            h = h01;
                        } else { // we are entirely in the interior or in the boundary
                            if (mesh_.valence(v0) < mesh_.valence(v1)) {
                                h = h01;
                            } else if (mesh_.valence(v0) > mesh_.valence(v1)) {
                                h = h10;
                            } else {
                                if (mesh_.is_collapse_ok(h01)) h = h01;
                                else if (mesh_.is_collapse_ok(h10)) h = h10;
                                else continue;
                            }
                        }

                        // collapse h if collapse is okay.
                        if (mesh_.is_collapse_ok(h)) {
                            finished = false;
                            mesh_.collapse(h);
                            c++;
                        }
                    }
                }
            }
        }

        cout << "Collapsed " << c << " short edges in " << i << " iterations." << endl;
        mesh_.garbage_collection();

        if (i == 100) std::cerr << "collapse break\n";
    }

    void MeshProcessing::equalize_valences() {
        Mesh::Edge_iterator e_it, e_end(mesh_.edges_end());
        Mesh::Vertex v0, v1, v2, v3;
        Mesh::Halfedge h;
        int ve0, ve1, ve2, ve3, ve_before, ve_after;
        bool finished;
        int i;
        int c = 0;

        // flip all edges
        for (finished = false, i = 0; !finished && i < 100; ++i) {
            finished = true;

            for (e_it = mesh_.edges_begin(); e_it != e_end; ++e_it) {
                if (!mesh_.is_boundary(*e_it) && !mesh_.is_deleted(*e_it)) {
                    // ------------- IMPLEMENT HERE ---------
                    //  Extract valences of the four vertices involved to an eventual flip.
                    //  Compute the sum of the squared valence deviances before flip
                    //  Compute the sum of the squared valence deviances after an eventual flip
                    //  If valence deviance is decreased and flip is possible, flip the vertex
                    //  Leave the loop running until no collapse has been done (use the finished variable)
                    // ------------- IMPLEMENT HERE ---------

                    // get the four vertices involved in a potential edge flip

                    auto h01 = mesh_.halfedge(*e_it, 0);
                    auto h10 = mesh_.opposite_halfedge(h01);

                    v0 = mesh_.to_vertex(h10);
                    v1 = mesh_.to_vertex(h01);
                    v2 = mesh_.to_vertex(mesh_.next_halfedge(h10));
                    v3 = mesh_.to_vertex(mesh_.next_halfedge(h01));

                    // Check if mesh has triangular structure and vertices are correct (only in Debug profile)
                    assert(mesh_.to_vertex(mesh_.next_halfedge(mesh_.next_halfedge(h10))) == v1);
                    assert(mesh_.to_vertex(mesh_.next_halfedge(mesh_.next_halfedge(h01))) == v0);

                    // compute deviations
                    auto compute_deviation = [this](Mesh::Vertex &v) -> int {
                        int opt = mesh_.is_boundary(v) ? 4 : 6;
                        return mesh_.valence(v) - opt;
                    };
                    ve0 = compute_deviation(v0);
                    ve1 = compute_deviation(v1);
                    ve2 = compute_deviation(v2);
                    ve3 = compute_deviation(v3);

                    // compute loss before and after
                    auto sqr = [](int i) -> int { return i * i; }; //convenience square function
                    ve_before = sqr(ve0) + sqr(ve1) + sqr(ve2) + sqr(ve3);
                    // v0 and v1 loose one edge and v2 and v3 gain one
                    // positive deviation means we have more edges than desired and vice versa.
                    ve_after = sqr(--ve0) + sqr(--ve1) + sqr(++ve2) + sqr(++ve3);

                    if (ve_after < ve_before && mesh_.is_flip_ok(*e_it)) {
                        finished = false;
                        mesh_.flip(*e_it);
                        c++;
                    }
                }
            }
        }

        cout << "Flipped " << c << " edges in " << i << " iterations." << endl;

        if (i == 100) std::cerr << "flip break\n";
    }

    void MeshProcessing::tangential_relaxation() {
        Mesh::Vertex_iterator v_it, v_end(mesh_.vertices_end());
        Point u, n;
        Point laplace;
        int c = 0;

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
                    laplace = Point(0.0);
                    for (const auto &v1 : mesh_.vertices(*v_it)) {
                        laplace += mesh_.position(v1) - mesh_.position(*v_it);
                    }
                    laplace /= mesh_.valence(*v_it);

                    // We use Eigen's geometry package to compute the projections
                    Eigen::Vector3d normal, laplacian, proj;
                    // convert surface_mesh vectors to Eigen vectors
                    n = normals[*v_it];
                    normal << n.x, n.y, n.z;
                    laplacian << laplace.x, laplace.y, laplace.z;

                    Eigen::Hyperplane<double, 3> tangentPlane(normal, 0);
                    proj = tangentPlane.projection(laplacian);

                    u = Point(proj(0), proj(1), proj(2));

                    // Set update to the projection
                    update[*v_it] = .1*u;
                }
            }

            // update the positions according to update
            for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it) {
                if (!mesh_.is_boundary(*v_it)) {
                    c++;
                    mesh_.position(*v_it) += update[*v_it];
                }
            }
        }
        cout << "Moved " << c << " vertices." << endl;
    }

    void MeshProcessing::calc_uniform_mean_curvature() {
        Mesh::Vertex_property <Scalar> v_unicurvature =
                mesh_.vertex_property<Scalar>("v:unicurvature", 0.0f);
        // ------------- IMPLEMENT HERE ---------
        // For each non-boundary vertex, approximate mean curvature using
        // the length of the uniform Laplacian approximation
        // Save your approximation in unicurvature vertex property of the mesh.
        // ------------- IMPLEMENT HERE ---------

        Point laplace(0.0);
        for (const auto &v : mesh_.vertices()) {
            laplace = 0;
            for (const auto &v2 : mesh_.vertices(v)) {
                laplace += mesh_.position(v2) - mesh_.position(v);
            }
            v_unicurvature[v] = norm(laplace / mesh_.valence(v));
        }
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

        Point laplace(0.0);
        for (const auto &v : mesh_.vertices()) {
            laplace = 0;
            for (const auto &v2 : mesh_.vertices(v)) {
                Mesh::Edge e = mesh_.find_edge(v, v2);
                laplace += e_weight[e] * (mesh_.position(v2) - mesh_.position(v));
            }
            laplace *= v_weight[v];
            v_curvature[v] = norm(laplace);
        }
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

        Mesh::Vertex_around_vertex_circulator vv_c, vv_c2;
        Scalar angles, cos_angle;
        Point d0, d1;
        Scalar lb(-.999f), ub(.999f);  // lower bound, upper bound for cos
        for (const auto &v : mesh_.vertices()) {
            if (mesh_.is_boundary(v)) continue;

            angles = 0.0f;

            vv_c = mesh_.vertices(v);
            vv_c2 = mesh_.vertices(v);

            surface_mesh::Vec3 pos = mesh_.position(v);
            for (const auto &v1 : vv_c) {
                ++vv_c2; //this is safe as the circulator iterator is implemented circularly. i.e. incrementing end -> begin

                d0 = mesh_.position(v1) - pos;
                d1 = mesh_.position(*vv_c2) - pos;

                cos_angle = min(ub, max(lb, dot(d0, d1) / (norm(d0) * norm(d1))));
                angles += acos(cos_angle);
            }

            v_gauss_curvature[v] = float(2 * M_PI - angles) * 2.0f * v_weight[v];
        }
    }

    void MeshProcessing::calc_max_curvature() {
        Mesh::Vertex_property <Scalar> v_gauss_curvature =
                mesh_.vertex_property<Scalar>("v:gauss_curvature", 0.0f);
        Mesh::Vertex_property <Scalar> v_curvature =
                mesh_.vertex_property<Scalar>("v:curvature", 0.0f);
        Mesh::Vertex_property <Scalar> v_max_curvature =
                mesh_.vertex_property<Scalar>("v:max_curvature", 0.0f);

        bool got_nans = false;
        for (const auto &v : mesh_.vertices()) {
            float H = v_curvature[v];
            float K = v_gauss_curvature[v];
            float k_max = H + sqrtf(H * H - K);
            if (k_max != k_max) { // we may get NaN --> use mean curvature as fallback
                got_nans = true;
                k_max = H;
            }
            v_max_curvature[v] = k_max;
        }
        if (got_nans) std::cerr << "Got NaNs computing max curvature -> using mean curvature for those" << endl;
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

            v_weight[v] = static_cast<float>(0.5 / area);
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
        Mesh::Vertex_property <Color> v_color_max_curv =
                mesh_.vertex_property<Color>("v:color_max_curv",
                                             Color(1.f, 1.f, 1.f));

        Mesh::Vertex_property <Color> v_color_target_length =
                mesh_.vertex_property<Color>("v:color_target_length",
                                             Color(1.f, 1.f, 1.f));

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
        Mesh::Vertex_property <Scalar> v_max_curvature =
                mesh_.vertex_property<Scalar>("v:max_curvature", 0.0f);

        calc_weights();
        calc_uniform_mean_curvature();
        calc_mean_curvature();
        calc_gauss_curvature();
        calc_max_curvature();
        color_coding(vertex_valence, &mesh_, v_color_valence, 3 /* min */,
                     8 /* max */);
        color_coding(v_unicurvature, &mesh_, v_color_unicurvature);
        color_coding(v_curvature, &mesh_, v_color_curvature);
        color_coding(v_gauss_curvature, &mesh_, v_color_gaussian_curv);

        color_coding(v_max_curvature, &mesh_, v_color_max_curv);

        // get the mesh attributes and upload them to the GPU
        int j = 0;
        unsigned int n_vertices(mesh_.n_vertices());

        // Create big matrices to send the data to the GPU with the required
        // format
        color_valence_ = Eigen::MatrixXf(3, n_vertices);
        color_unicurvature_ = Eigen::MatrixXf(3, n_vertices);
        color_curvature_ = Eigen::MatrixXf(3, n_vertices);
        color_gaussian_curv_ = Eigen::MatrixXf(3, n_vertices);
        color_max_curv_ = Eigen::MatrixXf(3, n_vertices);
        normals_ = Eigen::MatrixXf(3, n_vertices);
        points_ = Eigen::MatrixXf(3, n_vertices);
        indices_ = MatrixXu(3, mesh_.n_faces());
        color_target_length_ = Eigen::MatrixXf(3, n_vertices);

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

            color_max_curv_.col(j) << v_color_max_curv[v].x,
                    v_color_max_curv[v].y,
                    v_color_max_curv[v].z;
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


