#include <OpenGP/GL/TrackballWindow.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderShaded.h>
#include <OpenGP/SurfaceMesh/GL/SurfaceMeshRenderFlat.h>
#include <OpenGP/GL/PointsRenderer.h>
#include <OpenGP/GL/SegmentsRenderer.h>
#include <random>
#include <thread>

using namespace OpenGP;
using namespace std;
using namespace Eigen;

struct MainWindow : public TrackballWindow {
    PointsRenderer render_points = PointsRenderer();
    SegmentsRenderer render_segments = SegmentsRenderer();

    MatMxN points;
    MatMxN points_3d_render;
    int num_points;
    double radius = 0.3;
    SegmentsRenderer::Segments segments;

// ============================================================================
// Exercise 2 : fill the 2 functions below (see PDF for instructions)
// To test your implementation, use the S key for laplacian smoothing and the
// C key for the osculating circle.
// Hint : try to play with epsilon
// ============================================================================
    // time step for smoothing
    double epsilon = 0.001;

    double curve_length(const MatMxN &points_of_curve) {
        double curve_length = 0.;
        for (int i = 0; i < num_points; ++i) {
            curve_length += (points_of_curve.col(pyMod(i, num_points)) - points_of_curve.col(pyMod(i + 1, num_points))).norm();
        }
        return curve_length;
    }

    void rescale(MatMxN &points_updated) {
        points_updated *= curve_length(points) / curve_length(points_updated);
    }

    void laplacianSmoothing() {
        // Curve Smoothing - centroid (this function should do one iteration of smoothing)

        // create new blank matrix for updated points
        MatMxN points_updated;
        points_updated = MatMxN::Zero(2, num_points);

        // update second to second to last
        for (int i = 0; i < num_points; ++i) {
            points_updated.col(i) = (1. - epsilon) * points.col(i);
            points_updated.col(i) +=
                    epsilon * (points.col(pyMod(i - 1, num_points)) + points.col(pyMod(i + 1, num_points))) / 2.;
        }

        rescale(points_updated);
        points = points_updated;
    }

private:
    Vector2f computeCircleCenter(const Vector2f &x, const Vector2f &y, const Vector2f &z) {
        Matrix3f A, B, C;
        A << x(0), x(1), 1,
                y(0), y(1), 1,
                z(0), z(1), 1;
        B << x(0) * x(0) + x(1) * x(1), x(1), 1,
                y(0) * y(0) + y(1) * y(1), y(1), 1,
                z(0) * z(0) + z(1) * z(1), z(1), 1;
        C << x(0) * x(0) + x(1) * x(1), x(0), 1,
                y(0) * y(0) + y(1) * y(1), y(0), 1,
                z(0) * z(0) + z(1) * z(1), z(0), 1;
        Vector2f center(B.determinant() / (2 * A.determinant()),
                        -C.determinant() / (2 * A.determinant()));
        return center;
    }

    /**
     * In (plain) c++: -1 mod 6 = -1
     * In python: -1 mod 6 = 5
     *
     * As we want the python behaviour, we added this helper function
     */
    template<typename T, typename U>
    static T pyMod(const T &a, const U &b) {
        return (b + (a % b)) % b;
    }

public:
    void osculatingCircle() {
        // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
        MatMxN newPoints(points);
        for (int i = 0; i < points.cols(); i++) {
            Vector2f point = points.col(i);
            Vector2f delta = computeCircleCenter(points.col(pyMod(i - 1, points.cols())),
                                                 point,
                                                 points.col(pyMod(i + 1, points.cols())))
                             - point;
            newPoints.col(i) = point + epsilon * delta / delta.norm();
        }
        rescale(newPoints);
        points = newPoints;
    }

// ============================================================================
// END OF Exercise 2 (do not thouch the rest of the code)
// ============================================================================

    void generateRandomizedClosedPolyline() {
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0., 5 * 3e-2);

        Vec2 center(3e-2, 2e-3);


        points = MatMxN::Zero(2, num_points);
        for (int i = 0; i < num_points; ++i) {
            double frac = static_cast<double>(i) / static_cast<double>(num_points);
            points(0, i) = center(0) + radius * cos(2. * M_PI * frac) + distribution(generator);
            points(1, i) = center(1) + radius * sin(2. * M_PI * frac) + distribution(generator);
        }
    }

    void render() {

        // Prepare the render points
        points_3d_render = MatMxN::Zero(3, points.cols());
        points_3d_render.block(0, 0, 2, points.cols()) = points;

        // Rebuild the segments
        segments.clear();
        for (int i = 0; i < points_3d_render.cols(); ++i) {
            segments.push_back({points_3d_render.col(i), points_3d_render.col((i + 1) % points_3d_render.cols())});
        }
        render_points.init_data(points_3d_render);
        render_segments.init_data(segments);
    }

    MainWindow(int argc, char **argv) : TrackballWindow("2D Viewer", 640, 480) {
        num_points = 30;
        generateRandomizedClosedPolyline();

        this->scene.add(render_points);
        this->scene.add(render_segments);

        render();
    }

    bool key_callback(int key, int scancode, int action, int mods) override {
        TrackballWindow::key_callback(key, scancode, action, mods);
        if (key == GLFW_KEY_S && action == GLFW_RELEASE) {
            laplacianSmoothing();
        } else if (key == GLFW_KEY_C && action == GLFW_RELEASE) {
            osculatingCircle();
        } else if (key == GLFW_KEY_1 && action == GLFW_RELEASE) {
            num_points = 30;
            radius = 0.3;
            generateRandomizedClosedPolyline();
        } else if (key == GLFW_KEY_2 && action == GLFW_RELEASE) {
            num_points = 50;
            radius = 0.1;
            generateRandomizedClosedPolyline();
        } else if (key == GLFW_KEY_3 && action == GLFW_RELEASE) {
            num_points = 100;
            radius = 0.1;
            generateRandomizedClosedPolyline();
        } else if (key == GLFW_KEY_4 && action == GLFW_RELEASE) {
            num_points = 150;
            radius = 0.1;
            generateRandomizedClosedPolyline();
        }

        render();
        return true;
    }
};


int main(int argc, char **argv) {
    MainWindow window(argc, argv);
    return window.run();
}
