#include "../include/biharmonic_distance.h"
#include "../include/biharmonic_distance_approx.h"

#include <Eigen/Core>
#include <igl/isolines_map.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOFF.h>
#include <igl/unproject_onto_mesh.h>
#include <iostream>
#include <string>

using namespace Eigen;
using namespace std;

// This viusalization uses similar color map setting as tutorial #716
// for easy comparison between biharmonic and geodesic distances.
void set_colormap(igl::opengl::glfw::Viewer & viewer)
{
  // colormap texture
  int num_intervals = 50;
  MatrixXd CM(num_intervals, 3);
  for (int i = 0; i < num_intervals; i++) {
    double t = double(num_intervals - i - 1) / double(num_intervals - 1);
    CM(i, 0) = std::max(std::min(2.0 * t - 0.0, 1.0), 0.0);
    CM(i, 1) = std::max(std::min(2.0 * t - 1.0, 1.0), 0.0);
    CM(i, 2) = std::max(std::min(6.0 * t - 5.0, 1.0), 0.0);
  }
  igl::isolines_map(MatrixXd(CM), CM);
  viewer.data().set_colormap(CM);
}

int main(int argc, char *argv[])
{
  MatrixXd V;
  MatrixXi F;
  MatrixXd D;

  bool down_on_mesh = false;
  if ((argc != 1) && (argc != 2) && (argc != 3)) {
    cout<<"Usage:\n";
    cout<<"  exact computation: ./biharmonic_distance path_to_data\n";
    cout<<"  approximate computation: ./biharmonic_distance path_to_data k\n";
    cout<<"  (k is the number of eigenvectors to use for computing the approximated distance)\n";
    return 0;
  }

  if (argc == 2 || argc == 1) {
    cout<<"Use exact distance computation.\n";
    igl::read_triangle_mesh(
      (argc > 1 ? argv[1] : "../data/lucy-registration-complete.obj"), V, F);
    biharmonic_distance(V, F, D);
  }

  if (argc == 3) {
    char *p;
    int k = strtol(argv[2], & p, 10);
    cout<<"Use approximate distance computation.\n";
    igl::read_triangle_mesh(argv[1], V, F);
    biharmonic_distance_approx(V, F, k, D);
  }

  igl::opengl::glfw::Viewer viewer;
  const int xid = viewer.selected_data_index;
  viewer.append_mesh();

  // define the update function
  const auto update = [ & ]()->bool {
    int fid;
    Vector3f bc;
    // cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if (igl::unproject_onto_mesh(Vector2f(x, y), viewer.core().view,
        viewer.core().proj, viewer.core().viewport, V, F, fid, bc)) {
      VectorXd Dd;
      // if big mesh, just use closest vertex. Otherwise, blend distances to
      // vertices of face using barycentric coordinates.

      // 3d position of hit
      const RowVector3d m3 = V.row(F(fid, 0)) * bc(0) + V.row(F(fid, 1)) * bc(1) + V.row(F(fid, 2)) * bc(2);
      int cid = 0;
      Vector3d(
        (V.row(F(fid, 0)) - m3).squaredNorm(),
        (V.row(F(fid, 1)) - m3).squaredNorm(),
        (V.row(F(fid, 2)) - m3).squaredNorm()).minCoeff( & cid);
      const int vid = F(fid, cid);
      Dd = D.row(vid);

      viewer.data().set_data(Dd);
      return true;
    }
    return false;
  };

  // mouse call back function
  viewer.callback_mouse_down = [&](igl::opengl::glfw::Viewer & viewer, int, int)->bool {
    if (update()) {
      down_on_mesh = true;
      return true;
    }
    return false;
  };
  viewer.callback_mouse_move = [&](igl::opengl::glfw::Viewer & viewer, int, int)->bool {
    if (down_on_mesh) {
      update();
      return true;
    }
    return false;
  };
  viewer.callback_mouse_up = [&down_on_mesh](igl::opengl::glfw::Viewer & viewer, int, int)->bool {
    down_on_mesh = false;
    return false;
  };

  std::cout<<R"(
  drag mouse: change the source point of biharmonic distance.
  )";

  viewer.data().set_mesh(V, F);
  viewer.data().set_data(VectorXd::Zero(V.rows()));
  set_colormap(viewer);

  viewer.launch();
  update();
  return 0;
}
