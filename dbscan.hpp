#ifndef DBSCAN_HPP_
#define DBSCAN_HPP_

#include "kNeighbor.hpp"
#include <iostream>
#include <string>
#include <limits>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/io.h>
#include <glog/logging.h>

class Dbscan {
    typedef CGAL::Simple_cartesian<double> K;
	typedef K::Point_3 Point_3;
    typedef CGAL::Point_set_3<Point_3> Point_set;

public:
    Dbscan() {
        LOG(INFO) << "Init Dbscan";
        points_.clear();
        flag_.clear();
    };
    ~Dbscan() {};

    /**
    * @brief set input data
    * @param input_path format:*.asc 
    */
    bool set_input_points(std::string input_path);

    void set_eps(float eps) {eps_ = eps;};
    void set_minpts(int pts) {min_pts_ = pts;};

    bool compute();

    /**
     * @brief save files
     * 
     * @param minsize 
     * @return true 
     * @return false 
     */
    bool save_clusters(std::string path, int minsize);

    std::vector<std::vector<Point_3>> get_clusters() {
        return clusters_;
    };

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
    float eps_ = 15.0;
    int min_pts_ = 10;

    // expand cluster
    bool expand_cluster(const Point_3 &target, std::vector<Point_3> &cluster);
};

// @TODO add filesystem 
bool Dbscan::save_clusters(std::string path, int minsize = std::numeric_limits<int>::min()) {
    std::string output_file = "";
    for (int i = 0; i < clusters_.size(); ++i) {
        output_file = path + "/" + std::to_string(i) + "_.asc";
        std::cout << "output_file = " << output_file << std::endl;
        std::vector<Point_3> temp = clusters_[i];
        if (temp.size() > minsize) {
            std::ofstream out(output_file);
            for (int j = 0; j < temp.size(); j++)
                out << temp[j].x() << " " << temp[j].y() << " "  << temp[j].z() << std::endl;
            out.close();
        }        
    }
    return true;
}

bool Dbscan::expand_cluster(const Point_3 &target, std::vector<Point_3> &cluster) {
    ksearch_.searchK(Eigen::Vector3f(target.x(),target.y(),target.z()), min_pts_);
    std::vector<kNeighborData> knn = ksearch_.getQuary();
    if (knn.size() == 1) {
        flag_[knn[0].indices] = 0;
        return false;
    }
    if (knn[knn.size() - 1].square > eps_)
        return false;
    cluster.push_back(target);
    flag_[knn[0].indices] = 1;
    for (auto it : knn) {
        if (flag_[it.indices] == -1) {
            const Point_3 t(it.point(0), it.point(1), it.point(2));
            expand_cluster(t, cluster);
        }
    }
    return true;
}

bool Dbscan::compute() {
    // construct a kdtree
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
    std::vector<Eigen::Vector3f> data;
    for (int i = 0; i < points_.size(); ++i) {
        Eigen::Vector3f t = Eigen::Vector3f::Zero();
        t(0) = points_.point(i).x();
        t(1) = points_.point(i).y();
        t(2) = points_.point(i).z();
        data.push_back(t);
    } 
    ksearch_.setInputData(data);

    // add marked flag
    flag_.resize(points_.size(), -1);

    return true;
}

#endif // dbscan
