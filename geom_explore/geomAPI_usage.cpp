#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax3.hxx>
#include <Standard_Handle.hxx>
#include <Geom_Plane.hxx>
#include <Geom_SphericalSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <iostream>
#include <iomanip>

namespace geomAPI_examples {

  std::string formatValue(double value, double threshold = 1e-12, int precision = 6) {
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();  // Return the formatted value as a string
  }

  std::pair<double, double> projectPointOntoPlane(
    const gp_Pnt& point, 
    const gp_Ax3& planeAxis
  ) {
    // Define a plane using gp_Ax3 (origin, X-axis direction, and Z-axis direction)
    Handle(Geom_Plane) plane = new Geom_Plane(planeAxis);

    // Project the point onto the plane
    GeomAPI_ProjectPointOnSurf projector(point, plane);

    if (projector.IsDone() && projector.NbPoints() > 0) {
        // Get the projected point on the plane
        gp_Pnt projectedPoint = projector.Point(1);

        // Get (u, v) parameters of the projected point on the plane
        double u, v;
        projector.Parameters(1, u, v);

        // Output the results
        //std::cout << "Original Point: (" << point.X() << ", " << point.Y() << ", " << point.Z() << ")\n";
        //std::cout << "Projected Point: (" << projectedPoint.X() << ", " << projectedPoint.Y() << ", " << projectedPoint.Z() << ")\n";
        // std::cout << "\tU\t" << u << "\tV\t " << v;
        return std::pair<double, double>(u,v);
    } else {
        std::cerr << "Projection failed.\n";
        return std::make_pair(
          std::numeric_limits<double>::quiet_NaN(),
          std::numeric_limits<double>::quiet_NaN()
        );
    }
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

  gp_Pnt evalSphere(
    Handle(Geom_SphericalSurface) sphere,
    Standard_Real U, 
    Standard_Real V,
    bool printDerivs
  ) {
    // Evaluate the position and derivatives at the chosen (U, V)
    gp_Vec dU, dV, d2U, d2V, dUV;

    // Evaluate position and first derivatives (tangents)
    gp_Pnt position = sphere->Value(U, V);         // Point on the surface

    if(printDerivs) {
        std::cout << formatValue(U) << "\t" << formatValue(V) << "\t"
                  << formatValue(position.X()) << "\t"
                  << formatValue(position.Y()) << "\t" 
                  << formatValue(position.Z()) << std::endl;

        gp_Vec dU, dV;
        sphere->D1(U, V, position, dU, dV);    // First derivatives (tangents)
        std::cout << "\td/du = \t\t(" 
                  << formatValue(dU.X()) << ", " << formatValue(dU.Y()) << ", " << formatValue(dU.Z()) << ")" << std::endl;
        std::cout << "\td/dv = \t\t(" 
                  << formatValue(dV.X()) << ", " << formatValue(dV.Y()) << ", " << formatValue(dV.Z()) << ")" << std::endl;
        

        // Evaluate second derivatives (curvature)
        sphere->D2(U, V, position, dU, dV, d2U, d2V, dUV);  // Second derivatives
        std::cout << "\td2/du2 = \t(" 
                  << formatValue(d2U.X()) << ", " << formatValue(d2U.Y()) << ", " << formatValue(d2U.Z()) << ")" << std::endl;
        std::cout << "\td2/dv2 = \t(" 
                  << formatValue(d2V.X()) << ", " << formatValue(d2V.Y()) << ", " << formatValue(d2V.Z()) << ")" << std::endl;
        std::cout << "\td2/dudv = \t(" 
                  << formatValue(d2V.X()) << ", " << formatValue(dUV.Y()) << ", " << formatValue(dUV.Z()) << ")" << std::endl;
    }
    return position;
  };


  // this is an entry point into this file
  void sphere_example() {

    Handle(Geom_SphericalSurface) sphere;
    Standard_Real Umin, Umax, Vmin, Vmax;

    // Create an axis system (center and orientation of the sphere)
    gp_Ax3 axisSystem(gp_Pnt(0.0, 0.0, 0.0), gp_Dir(0.0, 0.0, 1.0));

    // Create a spherical surface with a radius of 10.0
    Standard_Real radius = 10.0;

    // leaks
    // const Geom_SphericalSurface* sphere = new Geom_SphericalSurface(axisSystem, radius);

    // does not leak
    //const std::unique_ptr<Geom_SphericalSurface> owner = std::make_unique<Geom_SphericalSurface>(axisSystem, radius);
    //const Geom_SphericalSurface* sphere = owner.get();

    // does not leak
    //const std::shared_ptr<Geom_SphericalSurface> owner = std::make_unique<Geom_SphericalSurface>(axisSystem, radius);
    //const Geom_SphericalSurface* sphere = owner.get();

    // does not leak
    sphere = new Geom_SphericalSurface(axisSystem, radius);

    /*
    Handle is a smart pointer used specifically within Open CASCADE for managing objects 
    derived from Standard_Transient and integrates with Open CASCADE’s memory management system, 
    while std::shared_ptr is a more general-purpose C++ smart pointer used across various types 
    of objects in standard C++ code.
    */

    // build a collection of surface positions
    const int MAX_I = 5000;
    std::array<gp_Pnt, MAX_I> pts; // collection of data to plot as points
    {
      // select some uvs on the surface
      std::array<std::pair<double, double>, MAX_I> uvs;
      for (int i = 0; i < MAX_I; i++) {
        uvs[i] = std::make_pair(
          0.0 +       5 * 2 * M_PI * i / MAX_I,  // winding around the sphere
          M_PI / 10 + 8 * 2 * M_PI * i / MAX_I   // like a LissaJous curve
        );
      }

      auto evaluator = [sphere](const std::pair<double, double>& uv) {
        return evalSphere(sphere, uv.first, uv.second, false);
      };
      // Use std::transform to apply the evaluation
      std::transform(
        uvs.begin(), 
        uvs.end(), 
        pts.begin(),
        evaluator
      );
    }

    // build a collection of surface positions
    const int MAX_LAT = 8;

    std::vector<std::array<gp_Pnt, MAX_I + 1>> linesVector; // collection of data to plot as polyline vxs
    for (int lineIndex = 0; lineIndex < MAX_LAT; lineIndex++) {
      std::array<gp_Pnt, MAX_I + 1> lines;
      std::array<std::pair<double, double>, MAX_I + 1> lines_uvs;
      for (int i = 0; i < MAX_I + 1; i++) {
        lines_uvs[i] = std::make_pair(
          2 * M_PI * i / MAX_I,   // winding around the sphere
          -M_PI / 2 + M_PI * lineIndex / MAX_LAT 
        );
      }
      auto evaluator = [sphere](const std::pair<double, double>& uv) {
        return evalSphere(sphere, uv.first, uv.second, false);
      };
      // Use std::transform to apply the evaluation
      std::transform(
        lines_uvs.begin(), 
        lines_uvs.end(), 
        lines.begin(),
        evaluator
      );
      linesVector.push_back(lines);
    }
    for (int lineIndex = 0; lineIndex < MAX_LAT; lineIndex++) {
      std::array<gp_Pnt, MAX_I + 1> lines;
      std::array<std::pair<double, double>, MAX_I + 1> lines_uvs;
      for (int i = 0; i < MAX_I + 1; i++) {
        lines_uvs[i] = std::make_pair(
          2 * M_PI * lineIndex / MAX_LAT,
          -M_PI / 2 + M_PI * i / MAX_I 
        );
      }
      auto evaluator = [sphere](const std::pair<double, double>& uv) {
        return evalSphere(sphere, uv.first, uv.second, false);
      };
      // Use std::transform to apply the evaluation
      std::transform(
        lines_uvs.begin(), 
        lines_uvs.end(), 
        lines.begin(),
        evaluator
      );
      linesVector.push_back(lines);
    }

    writeGeometryToJSON(pts, linesVector);
     
  }
}