#ifndef DBSCAN_HPP_
#define DBSCAN_HPP_

#include "kNeighbor.hpp"
#include <iostream>
#include <string>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/io.h>

class Dbscan {
    typedef CGAL::Simple_cartesian<double> K;
	typedef K::Point_3 Point_3;
    typedef CGAL::Point_set_3<Point_3> Point_set;
public:
    Dbscan();
    ~Dbscan();

    /**
    * @brief set input data
    * @param input_path format:*.asc 
    */
    bool set_input_points(std::string input_path);

    void set_eps(float eps) {eps_ = eps;};
    void set_minpts(int pts) {min_pts_ = pts;};

    bool compute();

private:
    Point_set points_;

    /**
     * @brief store each points attribute
     *         -1 : do not assigned
     *          0 : noisy
     *          1 : has assigned 
     */
    std::vector<int> flag_;

    std::vector<std::vector<Point_3>> clusters_;

    // class knn
    kNeighbor ksearch_;

    // eps && min_pts
    float eps_;
    int min_pts_;

    // expand cluster
    bool expand_cluster(const Point_3 &target, std::vector<Point_3> &cluster);
};

bool Dbscan::expand_cluster(const Point_3 &target, std::vector<Point_3> &cluster) {
    ksearch_.searchK(Eigen::Vector3f(target.x(),target.y(),target.z()), min_pts_);
    std::vector<Eigen::Vector4f> knn = ksearch_.getQuary();
    if (knn.size() == 1) {
        flag_[knn[0](3)] = 0;
        return false;
    }
    if (knn[knn.size() - 1](3) > eps_)
        return false; 
    cluster.push_back(target);
    flag_[knn[0](3)] = 1;
    for (auto it : knn) {
        if (flag_[it(3)] == -1) {
            const Point_3 t(it.x(), it.y(), it.z());
            expand_cluster(t, cluster);
        } 
    }
    return true;
}

bool Dbscan::compute() {
    for (int i = 0; i < points_.size(); ++i) {
        if (flag_[i] == -1) {
            std::vector<Point_3> cluster(0);
            expand_cluster(points_.point(i), cluster);
            clusters_.push_back(cluster);
        }
    }
    return true;
}

bool Dbscan::set_input_points(std::string input_path) {
    
    CGAL::IO::read_XYZ(input_path, points_);
    
    // add marked flag
    flag_.resize(points_.size(), -1);

    return true;
}

#endif // dbscan
