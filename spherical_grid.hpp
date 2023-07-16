#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
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

template <std::size_t CellCapacity = 10000U> struct SphericalCell
{
    static constexpr auto CELL_CAPACITY = CellCapacity;

    SphericalCell() : points{new SphericalPoint[CELL_CAPACITY]}, count(0U)
    {
    }

    ~SphericalCell()
    {
        delete[] points;
    }

    void reset()
    {
        count = 0U;
    }

    SphericalPoint *points;
    std::size_t count;
};

template <std::size_t CellCapacity = 10000U> class SphericalGrid
{
    static constexpr float TWO_PI = 2.0f * M_PIf;
    static constexpr float ONE_PI = M_PIf;
    static constexpr float HALF_PI = 0.5f * M_PIf;

  public:
    static constexpr auto CELL_CAPACITY = CellCapacity;

    explicit SphericalGrid(float radial_cell_size_m, float azimuthal_cell_size_rad, float elevation_cell_size_rad)
        : radial_cell_size_m_{radial_cell_size_m}, azimuthal_cell_size_rad_{azimuthal_cell_size_rad},
          elevation_cell_size_rad_{elevation_cell_size_rad}
    {
        radial_divisions_ = 1000U; // How to choose this?
        azimuthal_divisions_ = std::ceil(TWO_PI / azimuthal_cell_size_rad_);
        elevation_divisions_ = std::ceil(ONE_PI / elevation_cell_size_rad_);
        angular_divisions_ = azimuthal_divisions_ * elevation_divisions_;

        total_cells_ = radial_divisions_ * azimuthal_divisions_ * elevation_divisions_;
        spherical_grid_ = new SphericalCell<CELL_CAPACITY>[total_cells_];
    }

    ~SphericalGrid()
    {
        delete[] spherical_grid_;
    }

    void reset()
    {
        for (std::size_t i = 0U; i < total_cells_; ++i)
        {
            spherical_grid_->reset();
        }
    }

    inline std::size_t cellIndex(const SphericalPoint &point) noexcept
    {
        std::size_t radial_cell_index = static_cast<std::size_t>(point.radius_m / radial_cell_size_m_);

        // Wrap azimuthal angle to [0, 2*PI]
        float wrapped_azimuth = std::fmod(point.azimuth_rad, TWO_PI);
        if (wrapped_azimuth < 0.0f)
        {
            wrapped_azimuth += TWO_PI;
        }
        std::size_t azimuthal_cell_index = static_cast<std::size_t>(wrapped_azimuth / azimuthal_cell_size_rad_);

        // Wrap elevation angle to [0, PI]
        float wrapped_elevation = std::fmod(point.elevation_rad, ONE_PI);
        if (wrapped_elevation < 0.0f)
        {
            wrapped_elevation += ONE_PI;
        }
        std::size_t elevation_cell_index = static_cast<std::size_t>(wrapped_elevation / elevation_cell_size_rad_);

        std::size_t cell_index =
            radial_cell_index * angular_divisions_ + elevation_cell_index * azimuthal_divisions_ + azimuthal_cell_index;

        return cell_index;
    }

    inline void insert(const SphericalPoint &point)
    {
        std::size_t cell_index = cellIndex(point);
        if (spherical_grid_[cell_index].count >= CELL_CAPACITY)
        {
            throw std::runtime_error("Exceeded cell capacity");
        }
        spherical_grid_[cell_index].points[spherical_grid_[cell_index].count++] = point;
    }

    inline SphericalCell<CELL_CAPACITY> &getCell(const SphericalPoint &point)
    {
        std::size_t cell_index = cellIndex(point);
        return spherical_grid_[cell_index];
    }

    inline static float getDistance(const SphericalPoint &p1, const SphericalPoint &p2)
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
        const SphericalCell<CELL_CAPACITY> &cell = getCell(point);

        // Handle the case where the cell is empty
        bool found_neighbour = false;
        if (cell.count > 0)
        {
            nearest_neighbour.distance = getDistance(point, cell.points[0]);

            std::size_t nearest_neighbour_index = 0U;
            for (std::size_t i = 1; i < cell.count; ++i)
            {
                float distance = getDistance(point, cell.points[i]);
                if (distance < nearest_neighbour.distance)
                {
                    nearest_neighbour_index = i;
                    nearest_neighbour.distance = distance;
                }
            }

            nearest_neighbour.point = cell.points[nearest_neighbour_index];
            found_neighbour = true;
        }

        return found_neighbour;
    }

    struct CompareByDistance
    {
        inline bool operator()(const SphericalPoint &p, const SphericalPoint &p1, const SphericalPoint &p2)
        {
            return getDistance(p, p1) < getDistance(p, p2);
        }
    };

    inline bool getKNearestNeighbours(const SphericalPoint &point, std::size_t K,
                                      std::vector<SphericalPointWithDistance> &nearest_neighbours)
    {
        static constexpr CompareByDistance compare_by_distance{};

        nearest_neighbours.clear();

        std::size_t cell_index = cellIndex(point);
        std::vector<bool> visited(total_cells_, false);
        std::queue<std::size_t> to_visit;

        to_visit.push(cell_index);
        visited[cell_index] = true;

        while (!to_visit.empty())
        {
            std::size_t current_cell_index = to_visit.front();
            to_visit.pop();

            const SphericalCell<CELL_CAPACITY> &cell = spherical_grid_[current_cell_index];
            for (std::size_t i = 0; i < cell.count; ++i)
            {
                float d = getDistance(point, cell.points[i]);
                nearest_neighbours.push_back({cell.points[i], d});
            }

            if (nearest_neighbours.size() >= K)
            {
                std::partial_sort(nearest_neighbours.begin(), nearest_neighbours.begin() + K, nearest_neighbours.end(),
                                  compare_by_distance);
                nearest_neighbours.resize(K);
                return true;
            }

            for (std::size_t neighbor_index : cell.neighbors)
            {
                if (!visited[neighbor_index])
                {
                    to_visit.push(neighbor_index);
                    visited[neighbor_index] = true;
                }
            }
        }

        // If we reach this point, it means there are less than K points in the entire grid
        std::sort(nearest_neighbours.begin(), nearest_neighbours.end(), compare_by_distance);
        return false;
    }

  private:
    float radial_cell_size_m_;
    float azimuthal_cell_size_rad_;
    float elevation_cell_size_rad_;
    std::size_t radial_divisions_;
    std::size_t azimuthal_divisions_;
    std::size_t elevation_divisions_;
    std::size_t angular_divisions_;
    std::size_t total_cells_;
    SphericalCell<CELL_CAPACITY> *spherical_grid_;
};
} // namespace spatial_index