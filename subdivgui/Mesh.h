#include <igl/OpenGL_convenience.h>
#include <Eigen/Core>
#include <vector>
class Mesh
{
  public:
    // Rest and deformed positions
    Eigen::Matrix<double,Eigen::Dynamic,3> V,U,N,C;
    Eigen::MatrixXd W;
    Eigen::Matrix<float,Eigen::Dynamic,3,Eigen::RowMajor> UCT,NCT;
    // Quad and triangle face lists
    Eigen::Matrix<int,Eigen::Dynamic,4> Q;
    Eigen::Matrix<int,Eigen::Dynamic,3> F;
    Eigen::Matrix<int,Eigen::Dynamic,3,Eigen::RowMajor> T,T1;
    // Display list, vertex/normal/index buffer
    GLuint dl=0,vbo=0,nbo=0,ibo=0;
    // Display list, vbo, and nbo are stale
    bool wireframe = false;
    bool stale = true;
    bool must_update = true;
    bool force_visible;
    bool show_weights = false;
    std::vector<size_t> selected_weights;
    Mesh(bool _force_visible):force_visible(_force_visible){};
    void draw_and_cache();
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

