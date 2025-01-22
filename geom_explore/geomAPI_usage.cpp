#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax3.hxx>
#include <Standard_Handle.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <iostream>
#include <iomanip>

namespace geomAPI_examples {

  std::string formatValue(double value, double threshold = 1e-12, int precision = 6) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();  // Return the formatted value as a string
  }

  template <typename PtCollection, typename LineCollection>
  void writeGeometryToJSON(
    const PtCollection& pts,  // positions for points
    const LineCollection& polylines // positions for polyline vertices
  ) {
    try {
      // Write this data to a json file for viewing to pick up
      std::ofstream debugFile("viewer/public/output/view_data.json");
      if (!debugFile) {
          std::cerr << "Error opening file for writing.\n";
          return;
      }
      // we'll redirect cout to write to json file
      std::streambuf* coutBuf = std::cout.rdbuf(); // Save the original buffer

      std::cout.rdbuf(debugFile.rdbuf());  // Redirect std::cout to the file
      std::cout << "{\"pts\":[\n";
      for (size_t i = 0; i < pts.size(); ++i) {
          const gp_Pnt& pt = pts[i];
          std::cout 
            << "  { \"x\": " << pt.X() << ", \"y\": " << pt.Y() << ", \"z\": " << pt.Z() << "}";
          if (i != pts.size() - 1) {
            std::cout << ",";
          }
          std::cout << "\n"; // new point on a new line of the json file
      }
      std::cout << "],\n\"polylines\":[\n";

      for (size_t polyLineIndex = 0; polyLineIndex < polylines.size(); ++polyLineIndex) {
        std::cout << "[\n";
        const auto polyLine = polylines[polyLineIndex];
        for (size_t segIndex = 0; segIndex < polyLine.size(); ++segIndex) {

            const gp_Pnt& pt = polyLine[segIndex];
            std::cout 
              << "  { \"x\": " << pt.X() << ", \"y\": " << pt.Y() << ", \"z\": " << pt.Z() << "}";
            if (segIndex != polyLine.size() - 1) {
              std::cout << ",";
            }
            std::cout << "\n"; // new point on a new line of the json file
        }
        std::cout << "]";
        if (polyLineIndex != polylines.size() - 1) {
          std::cout << ",";
        }
        std::cout << "\n";
      }
      std::cout << "]}";
      // Restore std::cout to its original state
      std::cout.rdbuf(coutBuf);
    } catch (Standard_Failure& e) {
        std::cerr << "Error: " << e.GetMessageString() << std::endl;
    }  
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
    std::vector<gp_Pnt> pts(ptUvs.size());
    std::transform(ptUvs.begin(), ptUvs.end(), pts.begin(), evalSurface);

    std::vector<std::vector<gp_Pnt>> lineVxsVector; 
    for (const auto uvs : lineUVLists) {
      std::vector<gp_Pnt> lineVertices(uvs.size());
      std::transform(uvs.begin(), uvs.end(), lineVertices.begin(), evalSurface);
      lineVxsVector.push_back(lineVertices);
    }

    writeGeometryToJSON(pts, lineVxsVector);
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
      return sphere->Value(uv.first, uv.second);
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
    auto evalSphere = [torus](const std::pair<double, double>& uv) {
      return torus->Value(uv.first, uv.second);
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

    createOutputFile(ptUvs, lineUVLists, evalSphere);
  }

}