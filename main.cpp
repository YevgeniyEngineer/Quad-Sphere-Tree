#include "quad_sphere_tree.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <vector>

#define PI 3.14159265358979323846
#define NUM_POINTS 100000

using namespace spatial_index;

QuadSpherePoint randomPoint()
{
    return QuadSpherePoint{
        static_cast<double>(std::rand()) / RAND_MAX,              // radius_m in [0, 1]
        static_cast<double>(std::rand()) / RAND_MAX * 2 * PI,     // azimuth_rad in [0, 2*pi].
        static_cast<double>(std::rand()) / RAND_MAX * PI - PI / 2 // elevation_rad in [-pi/2, pi/2].
    };
}

bool operator==(const QuadSpherePoint &p1, const QuadSpherePoint &p2)
{
    static constexpr double tolerance = 1e-6;
    return (std::abs(p1.radius_m - p2.radius_m) < tolerance) &&
           (std::abs(p1.azimuth_rad - p2.azimuth_rad) < tolerance) &&
           (std::abs(p1.elevation_rad - p2.elevation_rad) < tolerance);
}

int main()
{
    // Seed the random number generator.
    std::srand(std::time(nullptr));

    // Create a tree.
    QuadSphereTree tree(NUM_POINTS);

    // Create a vector to hold the points for brute force search.
    std::vector<QuadSpherePoint> points;

    // Generate random points and insert them into the tree and the vector.
    for (std::size_t i = 0; i < NUM_POINTS; ++i)
    {
        QuadSpherePoint point = randomPoint();
        tree.insert(point);
        points.push_back(point);
    }

    // Generate a random target point.
    QuadSpherePoint target = randomPoint();

    // Search for the nearest neighbor using the tree.
    auto t1 = std::chrono::high_resolution_clock::now();

    QuadSpherePoint nearest_tree = tree.nearestNeighbourSearch(target);

    auto t2 = std::chrono::high_resolution_clock::now();

    // Search for the nearest neighbor using brute force.
    QuadSpherePoint nearest_brute;
    double min_distance = std::numeric_limits<double>::max();
    for (const QuadSpherePoint &point : points)
    {
        double distance = tree.distance(point, target);
        if (distance < min_distance)
        {
            min_distance = distance;
            nearest_brute = point;
        }
    }

    auto t3 = std::chrono::high_resolution_clock::now();

    // Check that the results are the same.
    if (nearest_tree == nearest_brute)
    {
        std::cout << "Test passed!" << std::endl;
    }
    else
    {
        std::cout << "Test failed!" << std::endl;
    }

    std::cout << "Elapsed time (tree search, seconds): " << static_cast<double>((t2 - t1).count()) * 1e-9 << std::endl;
    std::cout << "Elapsed time (brute force search, seconds): " << static_cast<double>((t3 - t2).count()) * 1e-9
              << std::endl;

    return 0;
}