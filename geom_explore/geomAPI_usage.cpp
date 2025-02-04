#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax3.hxx>
#include <Standard_Handle.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <iostream>
//#include <iomanip>
//#include <cmath>
#include "geom.hpp"

namespace geomAPI_examples {

  std::string formatValue(double value, double threshold = 1e-12, int precision = 6) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();  // Return the formatted value as a string
  }

  template <typename PointCollection, typename LineCollection, typename EvalFunc>
  void createOutputFile(
    const PointCollection& pt_uvs,
    const LineCollection& iso_v_lists,
    const LineCollection& iso_u_lists,
    const EvalFunc& eval_surface
  ) {
    std::vector<geom_examples::PtCollection> pts_colls;

    // build a collection of surface positions
    std::vector<geom_examples::Point> pts_to_display(pt_uvs.size());
    std::transform(pt_uvs.begin(), pt_uvs.end(), pts_to_display.begin(), eval_surface);

    pts_colls.push_back({
      .color = geom_examples::CYAN,
      .pts = pts_to_display,
      .isLine = false
    });

    for (const auto uvs : iso_v_lists) {
      std::vector<geom_examples::Point> line_vertices(uvs.size());
      std::transform(uvs.begin(), uvs.end(), line_vertices.begin(), eval_surface);
      pts_colls.push_back({
        .color = geom_examples::YELLOW,
        .pts = line_vertices,
        .isLine = true
      });
    }

    for (const auto uvs : iso_u_lists) {
      std::vector<geom_examples::Point> line_vertices(uvs.size());
      std::transform(uvs.begin(), uvs.end(), line_vertices.begin(), eval_surface);
      pts_colls.push_back({
        .color = geom_examples::RED,
        .pts = line_vertices,
        .isLine = true
      });
    }

    geom_examples::addGeometryToView(pts_colls);  
  }

  // this is an entry point into this file
  void sphereExample() {

    // Create an axis system for center and orientation of the sphere
    gp_Ax3 axis_system(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
    Standard_Real radius = 10.0;

    // leaks
    // const Geom_SphericalSurface* sphere = new Geom_SphericalSurface(axis_system, radius);

    // does not leak
    //const std::unique_ptr<Geom_SphericalSurface> owner = std::make_unique<Geom_SphericalSurface>(axis_system, radius);
    //const Geom_SphericalSurface* sphere = owner.get();

    // does not leak
    //const std::shared_ptr<Geom_SphericalSurface> owner = std::make_unique<Geom_SphericalSurface>(axis_system, radius);
    //const Geom_SphericalSurface* sphere = owner.get();

    // does not leak, uses internal OpenCascade memory management system
    Handle(Geom_SphericalSurface) sphere = new Geom_SphericalSurface(axis_system, radius);
    auto evalSphere = [sphere](const std::pair<double, double>& uv) {
      const gp_Pnt p = sphere->Value(uv.first, uv.second);
      return geom_examples::Point(p.X(), p.Y(), p.Z());
    };

    /*
    Handle is a smart pointer used specifically within Open CASCADE for managing objects 
    derived from Standard_Transient and integrates with Open CASCADE’s memory management system, 
    while std::shared_ptr is a more general-purpose C++ smart pointer used across various types 
    of objects in standard C++ code.
    */

    // select some uvs on the surface
    std::array<std::pair<double, double>, geom_examples::NUM_SAMPLES> pt_uvs;
    for (int i = 0; i < geom_examples::NUM_SAMPLES; i++) {
      pt_uvs[i] = std::make_pair(
        0.0 +       5 * 2 * M_PI * i / geom_examples::NUM_SAMPLES,  // winding around the sphere
        M_PI / 10 + 8 * 2 * M_PI * i / geom_examples::NUM_SAMPLES   // like a LissaJous curve
      );
    }

    std::vector<std::vector<std::pair<double, double>>> iso_v_lists;
    std::vector<std::vector<std::pair<double, double>>> iso_u_lists;

    // build uvs for isolines with constant V
    const int MAX_ISO = 8;
    for (int line_index = 0; line_index < MAX_ISO; line_index++) {
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, geom_examples::NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < geom_examples::NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * i / geom_examples::NUM_SAMPLES, // winding around the sphere
          -M_PI / 2 + M_PI * line_index / MAX_ISO // constant V
        );
      }
      iso_v_lists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }
    // build uvs for isolines with constant U
    for (int line_index = 0; line_index < MAX_ISO; line_index++) {
      std::array<gp_Pnt, geom_examples::NUM_SAMPLES + 1> line_vertices;
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, geom_examples::NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < geom_examples::NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * line_index / MAX_ISO, // constant U
          -M_PI / 2 + M_PI * i / geom_examples::NUM_SAMPLES // winding around the sphere
        );
      }
      iso_u_lists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }

    createOutputFile(pt_uvs, iso_v_lists, iso_u_lists, evalSphere);
  }

    // this is an entry point into this file
  void torusExample() {

    // Create an axis system for center and orientation of the sphere
    gp_Ax3 axis_system(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
    Standard_Real minor_radius = 5.0;
    Standard_Real major_radius = 30.0;

    Handle(Geom_ToroidalSurface) torus = new Geom_ToroidalSurface(axis_system, major_radius, minor_radius);
    auto evalTorus = [torus](const std::pair<double, double>& uv) {
      const gp_Pnt& p = torus->Value(uv.first, uv.second);
      return geom_examples::Point(p.X(), p.Y(), p.Z());
    };

    // select some uvs on the surface
    std::array<std::pair<double, double>, geom_examples::NUM_SAMPLES> pt_uvs;
    for (int i = 0; i < geom_examples::NUM_SAMPLES; i++) {
      pt_uvs[i] = std::make_pair(
        0.0 +       5 * 2 * M_PI * i / geom_examples::NUM_SAMPLES,  // winding around the sphere
        M_PI / 10 + 8 * 2 * M_PI * i / geom_examples::NUM_SAMPLES   // like a LissaJous curve
      );
    }

    std::vector<std::vector<std::pair<double, double>>> iso_v_lists;
    std::vector<std::vector<std::pair<double, double>>> iso_u_lists;

    // build uvs for isolines with constant V
    const int MAX_ISO = 8;
    for (int line_index = 0; line_index < MAX_ISO; line_index++) {
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, geom_examples::NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < geom_examples::NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * i / geom_examples::NUM_SAMPLES, // winding around the sphere
          2 * M_PI * line_index / MAX_ISO // constant V
        );
      }
      iso_v_lists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }
    // build uvs for isolines with constant U
    for (int line_index = 0; line_index < MAX_ISO; line_index++) {
      std::array<gp_Pnt, geom_examples::NUM_SAMPLES + 1> line_vertices;
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, geom_examples::NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < geom_examples::NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * line_index / MAX_ISO, // constant U
          2 * M_PI * i / geom_examples::NUM_SAMPLES // winding around the sphere
        );
      }
      iso_u_lists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }

    createOutputFile(pt_uvs, iso_v_lists, iso_u_lists, evalTorus);
  }

}