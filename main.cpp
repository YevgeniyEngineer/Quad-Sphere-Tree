#include "spherical_grid.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <limits>
#include <vector>

#define PI 3.14159265358979323846f
#define NUM_POINTS 150000

using namespace spatial_index;

SphericalPoint randomPoint()
{
    return SphericalPoint{
        static_cast<float>(std::rand()) / RAND_MAX,              // radius_m in [0, 1]
        static_cast<float>(std::rand()) / RAND_MAX * 2 * PI,     // azimuth_rad in [0, 2*pi].
        static_cast<float>(std::rand()) / RAND_MAX * PI - PI / 2 // elevation_rad in [-pi/2, pi/2].
    };
}

bool operator==(const SphericalPoint &p1, const SphericalPoint &p2)
{
    static constexpr float tolerance = 1e-5f;
    return (std::abs(p1.radius_m - p2.radius_m) < tolerance) &&
           (std::abs(p1.azimuth_rad - p2.azimuth_rad) < tolerance) &&
           (std::abs(p1.elevation_rad - p2.elevation_rad) < tolerance);
}

int main()
{
    // Seed the random number generator.
    std::srand(std::time(nullptr));

    // Create a grid.
    SphericalGrid grid(0.1f, 0.2f, 0.2f);

    // Create a vector to hold the points for brute force search.
    std::vector<SphericalPoint> points;
    for (std::size_t i = 0; i < NUM_POINTS; ++i)
    {
        SphericalPoint point = randomPoint();
        points.push_back(point);
    }

    // Generate random points and insert them into the grid and the vector.
    auto t00 = std::chrono::high_resolution_clock::now();
    for (const auto &point : points)
    {
        grid.insert(point);
    }
    auto t01 = std::chrono::high_resolution_clock::now();
    std::cout << "Elapsed time (grid insertion, seconds): " << static_cast<double>((t01 - t00).count()) * 1e-9
              << std::endl;

    // Generate a random target point.
    SphericalPoint target = randomPoint();

    // Search for the nearest neighbor using the grid.
    auto t1 = std::chrono::high_resolution_clock::now();
    SphericalPointWithDistance nearest_point;
    bool success = grid.getNearestNeighbor(target, nearest_point);
    auto t2 = std::chrono::high_resolution_clock::now();

    // Search for the nearest neighbor using brute force.
    SphericalPoint nearest_brute;
    float min_distance = std::numeric_limits<float>::max();
    for (const SphericalPoint &point : points)
    {
        float distance = grid.getDistance(point, target);
        if (distance < min_distance)
        {
            min_distance = distance;
            nearest_brute = point;
        }
    }
    auto t3 = std::chrono::high_resolution_clock::now();

    if (!success)
    {
        std::cout << "Could not find nearest neighbour using grid\n";
    }

    // Check that the results are the same.
    if (nearest_point.point == nearest_brute)
    {
        std::cout << "Test passed!" << std::endl;
    }
    else
    {
        std::cout << "Test failed!" << std::endl;
        std::cout << "Grid NN: radius_m = " << nearest_point.point.radius_m
                  << ", azimuth_rad = " << nearest_point.point.azimuth_rad
                  << ", elevation_rad = " << nearest_point.point.elevation_rad << std::endl;
        std::cout << "Brute-force NN: radius_m = " << nearest_brute.radius_m
                  << ", azimuth_rad = " << nearest_brute.azimuth_rad
                  << ", elevation_rad = " << nearest_brute.elevation_rad << std::endl;
    }

    std::cout << "Elapsed time (grid search, seconds): " << static_cast<double>((t2 - t1).count()) * 1e-9 << std::endl;
    std::cout << "Elapsed time (brute force search, seconds): " << static_cast<double>((t3 - t2).count()) * 1e-9
              << std::endl;

    return 0;
}
