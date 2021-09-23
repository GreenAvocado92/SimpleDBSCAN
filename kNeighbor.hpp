#pragma once
#ifndef KNEIGHBOR_HPP_
#define KNEIGHBOR_HPP_

// #define CGAL_LINKED_WITH_TBB 

#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/property_map.h>
#include <boost/iterator/zip_iterator.hpp>
#include <iostream>
#include <utility>
#include <vector>

#ifdef CGAL_LINKED_WITH_TBB 
	#include <tbb/blocked_range.h>
	#include <tbb/parallel_for.h>
#endif

using std::vector;
//glog
#include <glog/logging.h>

class kNeighbor
{
	typedef CGAL::Simple_cartesian<double> K;
	typedef K::Point_3 Point_d;
	typedef CGAL::Search_traits_3<K> TreeTraits;

	typedef boost::tuple<Point_d, int>  Point_and_int;
	typedef CGAL::Search_traits_adapter<Point_and_int,
		CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
		TreeTraits>                     Traits;
	typedef CGAL::Orthogonal_k_neighbor_search<Traits> Neighbor_search;
	typedef Neighbor_search::Tree Tree;

public:
	kNeighbor() {
		LOG(INFO) << "Init";
		points_.resize(0);
		indices_.resize(0);
		k_ = 1;
	};

	kNeighbor(std::vector<Point_d> points)
	{
		LOG(INFO) << "Init and setInputCloud";
		points_ = points;
		k_ = 1;

		build_tree_();
	}

	kNeighbor(std::vector<Eigen::Vector3f> points)
	{
		LOG(INFO) << "Init and setInputCloud";
		int i = 0;
		for (auto it : points) {
			points_.push_back(Point_d(it(0), it(1), it(2)));
			indices_.push_back(i);
			i++;
		}
		k_ = 1;

		build_tree_();
	}
	
	kNeighbor(std::vector<Eigen::Vector4f> points)
	{
		LOG(INFO) << "Init and setInputCloud";
		int i = 0;
		for (auto it : points) {
			points_.push_back(Point_d(it(0), it(1), it(2)));
			indices_.push_back(i);
			i++;
		}
		k_ = 1;

		build_tree_();
	}
	~kNeighbor() { };

	template <typename EigenT>
	void setInputData(std::vector<EigenT> points) {
		LOG(INFO) << "set input Data ";
		int i = 0;
		for (auto it : points) {
			points_.push_back(Point_d(it(0), it(1), it(2)));
			indices_.push_back(i);
			i++;
		}
		k_ = 1;

		build_tree_();
	}

	// search knn
	void searchK(Point_d t, size_t k);
	
	template <typename EigenT>
	void searchK(EigenT t, size_t k);
	
	// search rnn
	void searchR(Point_d t, float r);
	template <typename EigenT>
	void searchR(EigenT t, float r);

	// get quary 
	std::vector<Eigen::Vector4f> getQuary() { 
		LOG_IF(INFO, quary_.size() != 0) << "get knn quary successful";
		LOG_IF(INFO, quary_.size() == 0) << "get knn quary failed";
		return quary_; };
	
private:
	// 
	size_t k_;

	// store the raw point cloud
	std::vector<Point_d> points_;
	// the indices of points
	std::vector<int> indices_;

	// the result of knn or rnn
	std::vector<Eigen::Vector4f> quary_;

	// kdtree
	Tree tree_;

	bool isIndices_ = true;

	// init tree
	bool build_tree_();
};

template <typename EigenT>
void kNeighbor::searchK(EigenT t, size_t k)
{
	Point_d p(t(0), t(1), t(2));
	searchK(p, k);
}

void kNeighbor::searchK(Point_d t, size_t k)
{
	LOG(INFO) << "search K(" << k << ")nn  P(" << t << ")" ;
	quary_.clear();

	Neighbor_search search(tree_, t, k);

	for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
		int indices = boost::get<1>(it->first);
		Eigen::Vector4f point(boost::get<0>(it->first).x(), boost::get<0>(it->first).y(), boost::get<0>(it->first).z(), indices);

		quary_.push_back(point);
	}
}

void kNeighbor::searchR(Point_d t, float r)
{
	LOG(INFO) << "search R(" << r << ")nn  P(" << t << ")";
	quary_.clear();

	Neighbor_search search(tree_, t, 20);

	for (Neighbor_search::iterator it = search.begin(); it != search.end(); ++it) {
		int indices = boost::get<1>(it->first);
		Eigen::Vector4f point(boost::get<0>(it->first).x(), boost::get<0>(it->first).y(), boost::get<0>(it->first).z(), indices);

		if (it->second < r) quary_.push_back(point);
	}
}

bool kNeighbor::build_tree_()
{
	if (points_.size() == 0) {
		LOG(WARNING) << "No input points, tree build failed";
		return false;
	}
	if (!tree_.is_built()) {
		// if (isIndices_ == false)
		// 	tree_.insert(points_.begin(), points_.end());
		if (isIndices_ == true)
		 	tree_.insert(boost::make_zip_iterator(boost::make_tuple(points_.begin(), indices_.begin())),
		 		boost::make_zip_iterator(boost::make_tuple(points_.end(), indices_.end())));
		#ifdef CGAL_LINKED_WITH_TBB 
			tree_.build<CGAL::Parallel_tag>();
		#else
			tree_.build<CGAL::Sequential_tag>();
		#endif 
		LOG(INFO) << "tree build successful";
	}
	LOG(INFO) << "tree has built";
	return true;
}

template <typename EigenT>
void kNeighbor::searchR(EigenT t, float r)
{
	Point_d p(t(0), t(1), t(2));
	searchR(p, r);
}

#endif
