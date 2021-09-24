#include "dbscan.hpp"
#include <glog/logging.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) {
    google::InitGoogleLogging(argv[0]);
    FLAGS_logtostderr = false;
    FLAGS_alsologtostderr = false;
    FLAGS_stderrthreshold = google::ERROR;
 
    std::string input_path = argv[1];
    Dbscan db;
    db.set_input_points(input_path);
    db.compute();

    db.save_clusters(5000);

    return 0;
}