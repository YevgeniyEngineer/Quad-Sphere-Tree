#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <queue>
#include <stdexcept>
#include <vector>

namespace spatial_index
{
struct SphericalPoint
{
    float radius_m;
    float azimuth_rad;
    float elevation_rad;
};

struct SphericalPointWithDistance
{
    SphericalPoint point;
    float distance;
};

template <std::int32_t CellCapacity = 10000> struct SphericalCell
{
    static constexpr auto CELL_CAPACITY = CellCapacity;

    SphericalCell() : points{std::make_unique<SphericalPoint[]>(CELL_CAPACITY)}, count(0U)
    {
    }

    void reset()
    {
        count = 0U;
    }

    std::unique_ptr<SphericalPoint[]> points;
    std::int32_t count;
};

template <std::int32_t CellCapacity = 1000> class SphericalGrid
{
    static constexpr float TWO_PI = 2.0f * M_PIf;
    static constexpr float ONE_PI = M_PIf;
    static constexpr float HALF_PI = 0.5f * M_PIf;

  public:
    static constexpr auto CELL_CAPACITY = CellCapacity;

    static constexpr float MIN_RADIUS_M = 0.0f;
    static constexpr float MAX_RADIUS_M = 500.0f;
    static constexpr float MIN_AZIMUTH_RAD = 0.0f;
    static constexpr float MAX_AZIMUTH_RAD = TWO_PI;
    static constexpr float MIN_ELEVATION_RAD = 0.0f;
    static constexpr float MAX_ELEVATION_RAD = ONE_PI;

    explicit SphericalGrid(float radial_cell_size_m, float azimuthal_cell_size_rad, float elevation_cell_size_rad)
        : radial_cell_size_m_{radial_cell_size_m}, azimuthal_cell_size_rad_{azimuthal_cell_size_rad},
          elevation_cell_size_rad_{elevation_cell_size_rad}
    {
        radial_divisions_ = std::ceil((MAX_RADIUS_M - MIN_RADIUS_M) / radial_cell_size_m_);
        azimuthal_divisions_ = std::ceil((MAX_AZIMUTH_RAD - MIN_AZIMUTH_RAD) / azimuthal_cell_size_rad_);
        elevation_divisions_ = std::ceil((MAX_ELEVATION_RAD - MIN_ELEVATION_RAD) / elevation_cell_size_rad_);

        total_cells_ = radial_divisions_ * azimuthal_divisions_ * elevation_divisions_;
        std::cout << "Number of grid elements: " << total_cells_ << std::endl;

        spherical_grid_ = std::make_unique<SphericalCell<CELL_CAPACITY>[]>(total_cells_);
    }

    void reset()
    {
        for (std::int32_t i = 0U; i < total_cells_; ++i)
        {
            spherical_grid_[i].reset();
        }
    }

    struct CellIndices
    {
        std::int32_t radial;
        std::int32_t azimuthal;
        std::int32_t elevation;
    } cell_indices;

    void findCellIndices(const SphericalPoint &point) noexcept
    {
        // Set radius between [MIN_RADIUS_M, MAX_RADIUS_M]
        cell_indices.radial = std::min(
            radial_divisions_, static_cast<std::int32_t>((point.radius_m - MIN_RADIUS_M) / radial_cell_size_m_));

        // Wrap azimuthal angle to [MIN_AZIMUTH_RAD, MAX_AZIMUTH_RAD]
        static constexpr float total_azimuth_range_rad = MAX_AZIMUTH_RAD - MIN_AZIMUTH_RAD;
        float shifted_azimuth = point.azimuth_rad - MIN_AZIMUTH_RAD;
        float wrapped_azimuth = std::fmod(shifted_azimuth, total_azimuth_range_rad);
        if (wrapped_azimuth < 0.0f)
        {
            wrapped_azimuth += total_azimuth_range_rad;
        }
        cell_indices.azimuthal = static_cast<std::int32_t>(wrapped_azimuth / azimuthal_cell_size_rad_);

        // Wrap elevation angle to [MIN_ELEVATION_RAD, MAX_ELEVATION_RAD]
        static constexpr float total_elevation_range_rad = MAX_ELEVATION_RAD - MIN_ELEVATION_RAD;
        float shifted_elevation = point.elevation_rad - MIN_ELEVATION_RAD;
        float wrapped_elevation = std::fmod(shifted_elevation, total_elevation_range_rad);
        if (wrapped_elevation < 0.0f)
        {
            wrapped_elevation += total_elevation_range_rad;
        }
        cell_indices.elevation = static_cast<std::int32_t>(wrapped_elevation / elevation_cell_size_rad_);
    }

    inline std::int32_t cellIndex(std::int32_t radial_index, std::int32_t azimuthal_index, std::int32_t elevation_index)
    {
        return std::min((cell_indices.elevation * azimuthal_divisions_ + cell_indices.azimuthal) * radial_divisions_ +
                            cell_indices.radial,
                        total_cells_);
    }

    inline void insert(const SphericalPoint &point)
    {
        findCellIndices(point);

        auto &cell = spherical_grid_[cellIndex(cell_indices.radial, cell_indices.azimuthal, cell_indices.elevation)];
        cell.points[cell.count++] = point;
    }

    inline static float getDistanceSquared(const SphericalPoint &p1, const SphericalPoint &p2)
    {
        // d^2 = r1^2 + r2^2 - 2*r1*r2*(sin(ϕ1)*sin(ϕ2)*cos(θ1-θ2) + cos(ϕ1)*cos(ϕ2))

        float cos_angle =
            std::sin(p1.elevation_rad) * std::sin(p2.elevation_rad) * std::cos(p1.azimuth_rad - p2.azimuth_rad) +
            std::cos(p1.elevation_rad) * std::cos(p2.elevation_rad);

        float distance_squared =
            p1.radius_m * p1.radius_m + p2.radius_m * p2.radius_m - 2 * p1.radius_m * p2.radius_m * cos_angle;

        return distance_squared;
    }

    inline bool getNearestNeighbor(const SphericalPoint &point, SphericalPointWithDistance &nearest_neighbour)
    {
        findCellIndices(point);

        static constexpr std::int32_t RADIAL_SEARCH_OFFSET = 2;
        static constexpr std::int32_t AZIMUTHAL_SEARCH_OFFSET = 4;
        static constexpr std::int32_t ELEVATION_SEARCH_OFFSET = 4;

        std::int32_t min_radial_index = std::max(cell_indices.radial - RADIAL_SEARCH_OFFSET, 0);
        std::int32_t max_radial_index = std::min(radial_divisions_ - 1, cell_indices.radial + RADIAL_SEARCH_OFFSET);

        std::int32_t min_azimuthal_index = std::max(cell_indices.azimuthal - AZIMUTHAL_SEARCH_OFFSET, 0);
        std::int32_t max_azimuthal_index =
            std::min(azimuthal_divisions_ - 1, cell_indices.azimuthal + AZIMUTHAL_SEARCH_OFFSET);

        std::int32_t min_elevation_index = std::max(cell_indices.elevation - ELEVATION_SEARCH_OFFSET, 0);
        std::int32_t max_elevation_index =
            std::min(elevation_divisions_ - 1, cell_indices.elevation + ELEVATION_SEARCH_OFFSET);

        nearest_neighbour.distance = std::numeric_limits<float>::max();
        for (std::int32_t radial_index = min_radial_index; radial_index <= max_radial_index; ++radial_index)
        {
            for (std::int32_t azimuthal_index = min_azimuthal_index; azimuthal_index <= max_azimuthal_index;
                 ++azimuthal_index)
            {
                for (std::int32_t elevation_index = min_elevation_index; elevation_index <= max_elevation_index;
                     ++elevation_index)
                {
                    std::int32_t cell_index = cellIndex(radial_index, azimuthal_index, elevation_index);
                    const auto &cell = spherical_grid_[cell_index];

                    for (std::int32_t point_number = 0; point_number < cell.count; ++point_number)
                    {
                        const auto &cell_point = cell.points[point_number];
                        const float distance_squared = getDistanceSquared(point, cell_point);
                        if (distance_squared < nearest_neighbour.distance)
                        {
                            nearest_neighbour.distance = distance_squared;
                            nearest_neighbour.point = cell_point;
                        }
                    }
                }
            }
        }

        return true;
    }

    struct CompareByDistance
    {
        inline bool operator()(const SphericalPoint &p, const SphericalPoint &p1, const SphericalPoint &p2)
        {
            return getDistanceSquared(p, p1) < getDistanceSquared(p, p2);
        }
    };

  private:
    float radial_cell_size_m_;
    float azimuthal_cell_size_rad_;
    float elevation_cell_size_rad_;
    std::int32_t radial_divisions_;
    std::int32_t azimuthal_divisions_;
    std::int32_t elevation_divisions_;
    std::int32_t total_cells_;
    std::unique_ptr<SphericalCell<CELL_CAPACITY>[]> spherical_grid_;
};
} // namespace spatial_index