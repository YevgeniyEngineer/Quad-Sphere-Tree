#ifndef QUAD_SPHERE_TREE_HPP
#define QUAD_SPHERE_TREE_HPP

#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

namespace spatial_index
{
struct QuadSpherePoint
{
    double radius_m;
    double azimuth_rad;
    double elevation_rad;
};

struct QuadSphereNode
{
    QuadSpherePoint point;

    QuadSphereNode *north_west_ptr;
    QuadSphereNode *north_east_ptr;
    QuadSphereNode *south_west_ptr;
    QuadSphereNode *south_east_ptr;

    // Bounding box.
    // double min_radius_rad, max_radius_rad;
    // double min_azimuth_rad, max_azimuth_rad;
    // double min_elevation_rad, max_elevation_rad;

    void setPoint(const QuadSpherePoint &p)
    {
        point = p;
    }

    void setPoint(QuadSpherePoint &&p)
    {
        point = std::move(p);
    }

    QuadSphereNode()
        : north_west_ptr{nullptr}, north_east_ptr{nullptr}, south_west_ptr{nullptr}, south_east_ptr{nullptr}
    {
    }

    QuadSphereNode(const QuadSpherePoint &p)
        : point{p}, north_west_ptr{nullptr}, north_east_ptr{nullptr}, south_west_ptr{nullptr}, south_east_ptr{nullptr}
    {
    }
};

class QuadSphereTree
{
  public:
    QuadSphereTree() = delete;

    QuadSphereTree(std::size_t max_nodes)
        : max_nodes_{max_nodes}, node_buffer_{new QuadSphereNode[max_nodes_]}, next_free_node_index_{0U}, root_{nullptr}
    {
    }

    ~QuadSphereTree()
    {
        delete[] node_buffer_;
    }

    void insert(const QuadSpherePoint &point)
    {
        if (next_free_node_index_ == max_nodes_)
        {
            return;
        }

        insert(root_, point);
    }

    QuadSpherePoint nearestNeighbourSearch(const QuadSpherePoint &target_point)
    {
        QuadSpherePoint nearest_point;
        double min_distance = std::numeric_limits<double>::max();
        nearestNeighbourSearch(root_, target_point, nearest_point, min_distance);
        return nearest_point;
    }

    struct Haversine
    {
        inline double operator()(double value) const noexcept
        {
            return std::pow(std::sin(value / 2.0), 2.0);
        }
    };

    inline static double distance(const QuadSpherePoint &p1, const QuadSpherePoint &p2)
    {
        static constexpr Haversine haversine{};

        double d_azimuth_rad = p2.azimuth_rad - p1.azimuth_rad;
        double d_elevation_rad = p2.elevation_rad - p1.elevation_rad;

        double a = haversine(d_elevation_rad) + std::cos(p1.elevation_rad) * std::cos(p2.elevation_rad);
        double c = 2.0 * std::atan2(std::sqrt(a), std::sqrt(1.0 - a));

        return p1.radius_m * c;
    }

  private:
    std::size_t max_nodes_;
    QuadSphereNode *node_buffer_;
    std::size_t next_free_node_index_;
    QuadSphereNode *root_;

    static void insert(QuadSphereNode *&node, const QuadSpherePoint &point)
    {
        if (nullptr == node)
        {
            node = new QuadSphereNode(point);
            node->setPoint(point);
        }
        else
        {
            if (point.elevation_rad <= node->point.elevation_rad)
            {
                if (point.azimuth_rad <= node->point.azimuth_rad)
                {
                    insert(node->north_west_ptr, point);
                }
                else
                {
                    insert(node->north_east_ptr, point);
                }
            }
            else
            {
                if (point.azimuth_rad <= node->point.azimuth_rad)
                {
                    insert(node->south_west_ptr, point);
                }
                else
                {
                    insert(node->south_east_ptr, point);
                }
            }
        }
    }

    static void nearestNeighbourSearch(QuadSphereNode *node, const QuadSpherePoint &target_point,
                                       QuadSpherePoint &nearest_point, double &min_distance)
    {
        if (nullptr == node)
        {
            return;
        }

        // If the minimum possible distance to this node is greater than the best distance
        // found so far, we can skip searching this node.
        // double min_possible_distance = minPossibleDistance(target_point, node->min_azimuth_rad,
        // node->max_azimuth_rad,
        //                                                    node->min_elevation_rad, node->max_elevation_rad);
        // if (min_possible_distance >= min_distance)
        // {
        //     return;
        // }

        double current_distance = distance(node->point, target_point);

        if (current_distance < min_distance)
        {
            min_distance = current_distance;
            nearest_point = node->point;
        }

        nearestNeighbourSearch(node->north_west_ptr, target_point, nearest_point, min_distance);
        nearestNeighbourSearch(node->north_east_ptr, target_point, nearest_point, min_distance);
        nearestNeighbourSearch(node->south_west_ptr, target_point, nearest_point, min_distance);
        nearestNeighbourSearch(node->south_east_ptr, target_point, nearest_point, min_distance);
    }
};
} // namespace spatial_index

#endif // QUAD_SPHERE_TREE_HPP