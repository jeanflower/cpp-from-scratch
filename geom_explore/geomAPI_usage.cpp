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

  void sphere_example() {

    Handle(Geom_SphericalSurface) sphere;
    Standard_Real Umin, Umax, Vmin, Vmax;

    try {
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

      /*
      // Output some properties of the sphere
      std::cout << "Sphere radius: " << sphere->Radius() << std::endl;
      std::cout << "Sphere location: (" 
                << sphere->Position().Location().X() << ", " 
                << sphere->Position().Location().Y() << ", " 
                << sphere->Position().Location().Z() << ")" << std::endl;

      // Query the U and V ranges (parameterization)
      sphere->Bounds(Umin, Umax, Vmin, Vmax);

      // Print the U and V ranges
      std::cout << "U range: [" << Umin << ", " << Umax << "]" << std::endl;
      std::cout << "V range: [" << Vmin << ", " << Vmax << "]" << std::endl;
      */

    } catch (Standard_Failure& e) {
        std::cerr << "Error: " << e.GetMessageString() << std::endl;
    }

    /*
    try{
      // Choose a point (U, V) within the valid ranges
      Standard_Real U = Umin + (Umax - Umin) / 2;  // Choose middle of U range
      Standard_Real V = Vmin + (Vmax - Vmin) / 2;  // Choose middle of V range

      evalSphere(sphere, U, V, true);

    } catch (Standard_Failure& e) {
      std::cerr << "Error: " << e.GetMessageString() << std::endl;
    }
    */

    // Redirect std::cout to a file for viewing
    std::ofstream debugFile("viewer/public/output/xyz_coordinates.json");
    if (!debugFile) {
        std::cerr << "Error opening file for writing.\n";
        return;
    }

    // we'll redirect cout to write to json file
    std::streambuf* coutBuf = std::cout.rdbuf(); // Save the original buffer

    try{
      // select some uvs on the surface
      const int MAX_I = 3000;
      std::array<std::pair<double, double>, MAX_I> uvs;
      for (int i = 0; i < MAX_I; i++) {
        uvs[i] = std::make_pair(
          0.0 + 17 * 2 * M_PI * i / MAX_I, 
          0.1 + 11 * 2 * M_PI * i / MAX_I
        );
      }
      // build a collection of surface positions
      std::array<gp_Pnt, MAX_I> pts;
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

      std::cout.rdbuf(debugFile.rdbuf());          // Redirect std::cout to the file
      std::cout << "[\n";
      for (size_t i = 0; i < pts.size(); ++i) {
          const gp_Pnt& pt = pts[i];
          std::cout 
            << "  { \"x\": " 
            << pt.X() 
            << ", \"y\": " 
            << pt.Y() 
            << ", \"z\": " 
            << pt.Z() 
            << "}";
          if (i != pts.size() - 1) {
              std::cout << ",";
          }
          std::cout << "\n";
      }
      std::cout << "]\n";

    } catch (Standard_Failure& e) {
        std::cerr << "Error: " << e.GetMessageString() << std::endl;
    }

    // Restore std::cout to its original state
    std::cout.rdbuf(coutBuf);
  }
}