#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

namespace geom_examples {
  // Define a 3D point structure
  struct Point {
    double x, y, z;
    Point(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    double X() const {
      return x;
    }
    double Y() const {
      return y;
    }
    double Z() const {
      return z;
    }
  };

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
          const Point& pt = pts[i];
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

            const Point& pt = polyLine[segIndex];
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
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Unknown exception caught!" << std::endl;
    }
  }


  void indent(int level) {
    for (int i = 2; i > level; i--) {
      std::cout << "  ";
    }
  }

  // Define the NURBS curve class
  class NURBS {
    public:
      int degree;
      std::vector<Point> controlPoints;
      std::vector<double> weights;
      std::vector<double> knotVector;

      NURBS() : degree(0) {}

      NURBS(int degree, 
        const std::vector<Point>& controlPoints, 
        const std::vector<double>& weights, 
        const std::vector<double>& knotVector
      ): degree(degree), controlPoints(controlPoints), 
          weights(weights), knotVector(knotVector) {}

      // Cox-de Boor recursion formula to compute B-spline basis functions.
      // Include lots of comments as this is the most complex part of the code.
      // state as 'const' so that we can query the basisFunction of a 
      // const NURBS object
      double basisFunction(int i, int p, double t) const {
        if (p == 0) { // zero degree basis function
          // Base case: the 0th degree basis function is 1 if t is within the interval [t_i, t_{i+1}]
          const double result = (t >= knotVector[i] && t < knotVector[i + 1]) ? 1.0 : 0.0;
          //indent(p);
          //std::cout << "N(" << i << ", " << p << ", " << t << ") = " << result << std::endl;;
          return result;
        } else { // higher degree basis function
          // Recursive case
          //indent(p);
          //std::cout << "N(" << i << ", " << p << ", " << t << ") needs lower degree work ..." << std::endl;

          double denom1 = knotVector[i + p] - knotVector[i];
          double denom2 = knotVector[i + p + 1] - knotVector[i + 1];

          double num1 = (t - knotVector[i]);
          double num2 = (knotVector[i + p + 1] - t);

          double basis1 = basisFunction(i, p - 1, t);
          double basis2 = basisFunction(i + 1, p - 1, t);

          //indent(p);
          //std::cout << "N(" << i << ", " << p << ", " << t << ") = ";
          //std::cout << "(";

          double term1 = 0.0;
          if (denom1 > 1e-6) {
            term1 = num1 / denom1 * basis1;
            //std::cout << "(" << num1 << " / " << denom1 << " * " << basis1 << ")";
          } else {
            //std::cout << "0";
          }
          //std::cout << " + ";
          double term2 = 0.0;
          if (denom2 > 1e-6) {
            term2 = num2 / denom2 * basis2;
            //std::cout << "(" << num2 << " / " << denom2 << " * " << basis2 << ")";
          } else {
            //std::cout << "0";
          }

          //std::cout << ")";
          const double result = term1 + term2;
          //std::cout << " = " << result << std::endl << std::endl;
          return result;
        }
      }

      // Evaluate the NURBS curve at parameter t
      const Point evaluate(double t) const {
        int n = controlPoints.size();
        int p = degree;
        
        // Compute the weighted sum of control points
        double numeratorX = 0.0, numeratorY = 0.0, numeratorZ = 0.0;
        double denominator = 0.0;
        
        for (int i = 0; i <= n - 1; i++) {
          double N = basisFunction(i, p, t);
          // std::cout << "basisFunction("<<i<<", "<<p<<", "<<t<<") = " << N << std::endl;
          double weight = weights[i];
          //std::cout 
          //  << "parameter " << t 
          //  << " control point " << i 
          //  <<" basis function evaluation " << N 
          //  << " weight " << weight << std::endl << std::endl;
          numeratorX += N * weight * controlPoints[i].x;
          numeratorY += N * weight * controlPoints[i].y;
          numeratorZ += N * weight * controlPoints[i].z;
          denominator += N * weight;
        }

        // Normalize the result to get the actual curve point
        if (denominator != 0) {
          numeratorX /= denominator;
          numeratorY /= denominator;
          numeratorZ /= denominator;
        }

        return Point(numeratorX, numeratorY, numeratorZ);
      }
  };

  NURBS degree2_nurbs_example1() {
    NURBS nurbs = NURBS();
    {
      // Define degree, control points, weights, and knot vector
      int degree = 2;
      std::vector<Point> controlPoints = {
        Point(0, 0, 0), 
        Point(1, 1, 0), 
        Point(2, 0, 0)
      };
      std::vector<double> weights = {
        1.0, 
        1.0, 
        1.0
      };
      std::vector<double> knotVector = {
        0.0,
        0.0,
        0.0,
        // 0.5, too many knots! Fails reading control points out-of-range
        1.0,
        1.0,
        1.0
      };
      // Create the NURBS curve
      nurbs = NURBS(degree, controlPoints, weights, knotVector);
    }
    return nurbs;
  }

  void nurbs_example() {

    /*
    // check for read out of bounds
    std::vector<double> weights = {
      1.0, 
      1.0, 
      1.0
    };
    std::cout << "weights[3] = " << weights[3] << std::endl;
    */
   const NURBS nurbs = degree2_nurbs_example1();

    /*
    // Basis functions should sum to 1
    const double t = 0.499;
    double sum = 0.0;
    for (int i = 0; i < nurbs.controlPoints.size(); ++i) {
      const double N = nurbs.basisFunction(i, 2, t);
      std::cout << "N = " << N << std::endl;
      sum += N;
    }
    std::cout << "Sum of basis functions at t = " << t << ": " << sum << std::endl;
    */

    // Evaluate the NURBS curve at a chosen parameter value
    double t = 0.5;
    Point result = nurbs.evaluate(t);
    std::cout << "NURBS Curve evaluated at t = " << t << ": (" 
              << result.x << ", " << result.y << ", " << result.z << ")" << std::endl;

    // evaluate at a set of points and display
    std::vector<Point> ptsToDisplay;
    const int NUM_SAMPLES = 30;
    const double start = nurbs.knotVector[0];
    const double end = nurbs.knotVector[nurbs.knotVector.size() - 1];
    for (int i = 0; i < NUM_SAMPLES; i++) {
      const double t = start + i * (end - start) / (NUM_SAMPLES - 1);
      const Point& result = nurbs.evaluate(t);
      ptsToDisplay.push_back(Point(result.x, result.y, result.z));
      // std::cout << t << ", " 
      // << result.x << ", " << result.y << ", " << result.z << std::endl;
    }
    std::vector<std::vector<Point>> linesToDisplay;
    writeGeometryToJSON(ptsToDisplay, linesToDisplay);
  }

}