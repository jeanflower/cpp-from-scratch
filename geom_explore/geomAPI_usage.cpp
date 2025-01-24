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


  // instead of spelling out (and fixing) the types here, use a template
  // and let the compiler work out how to map the types
  template <typename PtCollection, typename LineCollection, typename EvalFunc>
  void createOutputFile(
    const PtCollection& ptUvs,
    const LineCollection& lineUVLists,
    const EvalFunc& evalSurface
  ) {
    // build a collection of surface positions
    std::vector<geom_examples::Point> pts(ptUvs.size());
    std::transform(ptUvs.begin(), ptUvs.end(), pts.begin(), evalSurface);

    std::vector<std::vector<geom_examples::Point>> lineVxsVector; 
    for (const auto uvs : lineUVLists) {
      std::vector<geom_examples::Point> lineVertices(uvs.size());
      std::transform(uvs.begin(), uvs.end(), lineVertices.begin(), evalSurface);
      lineVxsVector.push_back(lineVertices);
    }

    geom_examples::writeGeometryToJSON(pts, lineVxsVector);
  }

  // this is an entry point into this file
  void sphere_example() {

    // Create an axis system for center and orientation of the sphere
    gp_Ax3 axisSystem(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
    Standard_Real radius = 10.0;

    // leaks
    // const Geom_SphericalSurface* sphere = new Geom_SphericalSurface(axisSystem, radius);

    // does not leak
    //const std::unique_ptr<Geom_SphericalSurface> owner = std::make_unique<Geom_SphericalSurface>(axisSystem, radius);
    //const Geom_SphericalSurface* sphere = owner.get();

    // does not leak
    //const std::shared_ptr<Geom_SphericalSurface> owner = std::make_unique<Geom_SphericalSurface>(axisSystem, radius);
    //const Geom_SphericalSurface* sphere = owner.get();

    // does not leak, uses internal OpenCascade memory management system
    Handle(Geom_SphericalSurface) sphere = new Geom_SphericalSurface(axisSystem, radius);
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

    const int NUM_SAMPLES = 5000;

    // select some uvs on the surface
    std::array<std::pair<double, double>, NUM_SAMPLES> ptUvs;
    for (int i = 0; i < NUM_SAMPLES; i++) {
      ptUvs[i] = std::make_pair(
        0.0 +       5 * 2 * M_PI * i / NUM_SAMPLES,  // winding around the sphere
        M_PI / 10 + 8 * 2 * M_PI * i / NUM_SAMPLES   // like a LissaJous curve
      );
    }

    std::vector<std::vector<std::pair<double, double>>> lineUVLists;

    // build uvs for isolines with constant V
    const int MAX_ISO = 8;
    for (int lineIndex = 0; lineIndex < MAX_ISO; lineIndex++) {
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * i / NUM_SAMPLES, // winding around the sphere
          -M_PI / 2 + M_PI * lineIndex / MAX_ISO // constant V
        );
      }
      lineUVLists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }
    // build uvs for isolines with constant U
    for (int lineIndex = 0; lineIndex < MAX_ISO; lineIndex++) {
      std::array<gp_Pnt, NUM_SAMPLES + 1> lineVertices;
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * lineIndex / MAX_ISO, // constant U
          -M_PI / 2 + M_PI * i / NUM_SAMPLES // winding around the sphere
        );
      }
      lineUVLists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }

    createOutputFile(ptUvs, lineUVLists, evalSphere);
  }

    // this is an entry point into this file
  void torus_example() {

    // Create an axis system for center and orientation of the sphere
    gp_Ax3 axisSystem(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));
    Standard_Real minorRadius = 5.0;
    Standard_Real majorRadius = 30.0;

    Handle(Geom_ToroidalSurface) torus = new Geom_ToroidalSurface(axisSystem, majorRadius, minorRadius);
    auto evalTorus = [torus](const std::pair<double, double>& uv) {
      const gp_Pnt& p = torus->Value(uv.first, uv.second);
      return geom_examples::Point(p.X(), p.Y(), p.Z());
    };

    const int NUM_SAMPLES = 5000;

    // select some uvs on the surface
    std::array<std::pair<double, double>, NUM_SAMPLES> ptUvs;
    for (int i = 0; i < NUM_SAMPLES; i++) {
      ptUvs[i] = std::make_pair(
        0.0 +       5 * 2 * M_PI * i / NUM_SAMPLES,  // winding around the sphere
        M_PI / 10 + 8 * 2 * M_PI * i / NUM_SAMPLES   // like a LissaJous curve
      );
    }

    std::vector<std::vector<std::pair<double, double>>> lineUVLists;

    // build uvs for isolines with constant V
    const int MAX_ISO = 8;
    for (int lineIndex = 0; lineIndex < MAX_ISO; lineIndex++) {
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * i / NUM_SAMPLES, // winding around the sphere
          2 * M_PI * lineIndex / MAX_ISO // constant V
        );
      }
      lineUVLists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }
    // build uvs for isolines with constant U
    for (int lineIndex = 0; lineIndex < MAX_ISO; lineIndex++) {
      std::array<gp_Pnt, NUM_SAMPLES + 1> lineVertices;
      // closed isolines have repeated start/end uvs, so NUM_SAMPLES + 1
      std::array<std::pair<double, double>, NUM_SAMPLES + 1> uvs;
      for (int i = 0; i < NUM_SAMPLES + 1; i++) {
        uvs[i] = std::make_pair(
          2 * M_PI * lineIndex / MAX_ISO, // constant U
          2 * M_PI * i / NUM_SAMPLES // winding around the sphere
        );
      }
      lineUVLists.push_back(std::vector<std::pair<double, double>>(uvs.begin(), uvs.end()));
    }

    createOutputFile(ptUvs, lineUVLists, evalTorus);
  }

}