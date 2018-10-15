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
    double epsilon = 0.01;

    double curve_length(MatMxN points_of_curve) {
        double curve_length = 0.;

        for (int i = 0; i < num_points - 1; ++i) {
            curve_length += (points.col(i) - points.col(i+1)).norm();
        }
        curve_length += (points.col(num_points - 1) - points.col(0)).norm();

        cout << curve_length << endl;

        return curve_length;
    }

    void rescale(MatMxN &points_updated) {
        double original_length = curve_length(points);
        double temp_length = curve_length(points_updated);

        double ratio = original_length / temp_length;

        points_updated *= ratio;
    }

    void laplacianSmoothing() {
        // Curve Smoothing - centroid (this function should do one iteration of smoothing)

        // create new blank matrix for updated points
        MatMxN points_updated;
        points_updated = MatMxN::Zero(2, num_points);

        // update first
        points_updated.col(0) = (1.-epsilon) * points.col(0);
        points_updated.col(0) += epsilon * (points.col(num_points - 1) + points.col(1)) / 2.;

        // update second to second to last
        for (int i = 1; i < num_points - 1; ++i) {
            points_updated.col(i) = (1.-epsilon) * points.col(i);
            points_updated.col(i) += epsilon * (points.col(i-1) + points.col(i+1)) / 2.;
        }

        // update last
        points_updated.col(num_points - 1) = (1.-epsilon) * points.col(num_points - 1);
        points_updated.col(num_points - 1) += epsilon * (points.col(num_points - 2) + points.col(0)) / 2.;

        rescale(points_updated);
        points = points_updated;
    }

private:
    Vector2f computeCircleCenter(const Vector2f &x, const Vector2f &y, const Vector2f &z) {
        float ma = (y(1) - x(1)) / (y(0) - x(0));
        float mb = (z(1) - z(1)) / (z(0) - y(0));
        float cx = (ma * mb * (x(1) - z(1)) + mb * (x(0) - y(0)) - ma * (y(0) + z(0))) / (2 * (mb - ma));
        float cy = (-1 / ma) * (cx - (x(0) + z(0)) * 0.5) + (x(1) + y(1)) * 0.5;
        return Vector2f(cx, cy);
    }

    template<typename T, typename U>
    T pythonMod(const T &a, const U &b) {
        return (b + (a % b)) % b;
    }

public:
    void osculatingCircle() {
        // Curve Smoothing - osculating circle (again, this function should do one iteration of smoothing)
        MatMxN newPoints(points);
        for (int i = 0; i < points.cols(); i++) {
            Vector2f point = points.col(i);
            Vector2f delta = computeCircleCenter(points.col(pythonMod(i - 1, points.cols())),
                                                 point,
                                                 points.col(pythonMod(i + 1, points.cols())))
                             - point;
            delta.normalize();
            newPoints.col(i) = point + epsilon * delta;
        }
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

    std::shared_ptr<std::thread> t;

    MainWindow(int argc, char **argv) : TrackballWindow("2D Viewer", 640, 480) {
        num_points = 30;
        generateRandomizedClosedPolyline();

        this->scene.add(render_points);
        this->scene.add(render_segments);

        render();

        t = std::make_shared<std::thread>([this]() {
            while (false) {
                sleep(1);
                this->osculatingCircle();
                this->render();
                cout << "Update" << endl;
            }
        });
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
