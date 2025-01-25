#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <tuple>
#include <map>

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
  double basisFunctionNoCache(const std::vector<double>& knotVector, int i, int p, double t) {
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

      double basis1 = basisFunctionNoCache(knotVector, i, p - 1, t);
      double basis2 = basisFunctionNoCache(knotVector, i + 1, p - 1, t);

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

  struct KeyHash {
    std::size_t operator()(const std::tuple<int, int, double>& key) const {
      auto [i, p, t] = key;
      std::size_t h1 = std::hash<int>{}(i);
      std::size_t h2 = std::hash<int>{}(p);
      std::size_t h3 = std::hash<double>{}(t);
      
      // Use a more robust hash combination
      std::size_t hash = h1;
      hash ^= h2 + 0x9e3779b9 + (hash << 6) + (hash >> 2);  // Combine h1 and h2
      hash ^= h3 + 0x9e3779b9 + (hash << 6) + (hash >> 2);  // Combine hash and h3
      return hash;
    }
  };

  class BasisEvaluator {
    private:
      std::vector<double> knotVector;

    public:
      int cacheStrategy = 0;
      // strategy 0 uses no cache
      std::map<std::tuple<int, int, double>, double> basisCache1; // Cache for strategy 1
      std::unordered_map<std::tuple<int, int, double>, double, KeyHash> basisCache2; // Cache for strategy 2

      BasisEvaluator(
        const std::vector<double>& knotVector,
        int cacheStrategy = 0
      ): knotVector(knotVector), cacheStrategy(cacheStrategy) {}

      double basisFunction(int i, int p, double t) {
        if (cacheStrategy == 0) {
          return basisFunctionNoCache(knotVector, i, p, t);
        } else if(cacheStrategy == 1) {
          // Use a std::map as a cache to store the basis function values
          auto key = std::make_tuple(i, p, t);

          // Use emplace to avoid redundant lookups
          auto [iter, inserted] = basisCache1.emplace(key, 0.0); // Insert a default value (0.0) if key is missing
          if (inserted) {
              // Cache miss, compute the value
              iter->second = basisFunctionNoCache(knotVector, i, p, t);
              // Optionally log the cache miss
              // std::cout << "Cache miss for basis function (" << i << ", " << p << ", " << t << ")" << std::endl;
          } else {
              // Optionally log the cache hit
              // std::cout << "Cache hit for basis function (" << i << ", " << p << ", " << t << ")" << std::endl;
          }
          // Return the cached value
          return iter->second;
        } else if (cacheStrategy == 2) {
            // Use an unordered_map-based cache for faster lookups
            auto key = std::make_tuple(i, p, t);

            // Use emplace to avoid redundant lookups
            auto [iter, inserted] = basisCache2.emplace(key, 0.0); // Insert a default value (0.0) if key is missing
            if (inserted) {
                // Cache miss, compute the value
                iter->second = basisFunctionNoCache(knotVector, i, p, t);
            }
            // Return the cached value
            return iter->second;
          } else {
          return basisFunctionNoCache(knotVector, i, p, t);
        }
      }
  };

  // Define the NURBS curve class
  class NURBS {
    public:
      int degree;
      std::vector<Point> controlPoints;
      std::vector<double> weights;
      std::vector<double> knotVector;

    private:
      bool useCache = false;
      mutable BasisEvaluator basisEvaluator;

    public:
      NURBS(
        int degree, 
        const std::vector<Point>& controlPoints, 
        const std::vector<double>& weights, 
        const std::vector<double>& knotVector,
        int cacheStrategy
      ): degree(degree), controlPoints(controlPoints), 
          weights(weights), knotVector(knotVector), 
          basisEvaluator(knotVector, cacheStrategy) {
        if (degree + controlPoints.size() + 1 != knotVector.size()) {
          std::cerr << "Invalid NURBS curve definition: degree + control points + 1 != knot vector size" << std::endl;
        }
        if (weights.size() != controlPoints.size()) {
          std::cerr << "Invalid NURBS curve definition: weights size != control points size" << std::endl;
        }
      }

      double basisFunction(int i, int p, double t) const {
        return basisEvaluator.basisFunction(i, p, t);
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

  NURBS degree2_nurbs_example1(
    int cacheStrategy = 0
  ) {
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
    std::vector<double> knotVector = {  // 3 + 2 + 1 = 6 knots
      0.0,
      0.0,
      0.0,
      // 0.5, too many knots! Fails reading control points out-of-range
      1.0,
      1.0,
      1.0
    };
    BasisEvaluator basisEvaluator(knotVector, cacheStrategy);
    // Create the NURBS curve
    NURBS nurbs(degree, controlPoints, weights, knotVector, cacheStrategy);
    return nurbs;
  }

  NURBS degree3_nurbs_example1(
    int cacheStrategy = 0
  ) {
    // Define degree, control points, weights, and knot vector
    int degree = 3;
    std::vector<Point> controlPoints = {
      Point(0, 0, 0), 
      Point(1, 1, 0), 
      Point(2, 0, 0),
      Point(3, 1, 0), 
    };
    std::vector<double> weights = {
      1.0, 
      1.0, 
      1.0,
      1.0,
    };
    std::vector<double> knotVector = {   // 4 + 3 + 1 = 8 knots
      0.0,
      0.0,
      0.0,
      0.2,
      0.8,
      1.0,
      1.0,
      1.0
    };
    // Create the NURBS curve
    NURBS nurbs(degree, controlPoints, weights, knotVector, cacheStrategy);
    return nurbs;
  }

  NURBS degree4_nurbs_example1(
    int cacheStrategy = 0
  ) {
    // Define degree, control points, weights, and knot vector
    int degree = 4;
    std::vector<Point> controlPoints = {
      Point(0, 0, 0), 
      Point(1, 1, 0), 
      Point(2, 0, 0),
      Point(3, 1, 0), 
      Point(4, 0, 0),
    };
    std::vector<double> weights = {
      1.0, 
      1.0, 
      1.0,
      1.0,
      1.0,
    };
    std::vector<double> knotVector = { // 5 + 4 + 1 = 10 knots
      0.0,
      0.0,
      0.0,
      0.2,
      0.4,
      0.6,
      0.8,
      1.0,
      1.0,
      1.0
    };
    // Create the NURBS curve
    NURBS nurbs(degree, controlPoints, weights, knotVector, cacheStrategy);
    return nurbs;
  }

  NURBS degree6_nurbs_example1(
    int cacheStrategy = 0
  ) {
    // Define degree, control points, weights, and knot vector
    int degree = 6;
    std::vector<Point> controlPoints = {
      Point(0, 0, 0), 
      Point(1, 1, 0), 
      Point(2, 0, 0),
      Point(3, 1, 0), 
      Point(4, 0, 0),
      Point(5, 1, 0), 
      Point(6, 0, 0),
    };
    std::vector<double> weights = {
      1.0, 
      1.0, 
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
    };
    std::vector<double> knotVector = { // 7 + 6 + 1 = 14 knots
      0.0,
      0.0,
      0.0,
      0.2,
      0.2,
      0.4,
      0.4,
      0.6,
      0.6,
      0.8,
      0.8,
      1.0,
      1.0,
      1.0
    };
    // Create the NURBS curve
    NURBS nurbs(degree, controlPoints, weights, knotVector, cacheStrategy);
    return nurbs;
  }

  void nurbs_example(
    bool useCache = false
  ) {

    /*
    // check for read out of bounds
    std::vector<double> weights = {
      1.0, 
      1.0, 
      1.0
    };
    std::cout << "weights[3] = " << weights[3] << std::endl;
    */
   const NURBS nurbs = degree2_nurbs_example1(useCache);

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

  double doConfiguredWork(
    const NURBS& nurbs,
    const int numSamples, 
    const int numEvals
  ) {
    auto start = std::chrono::high_resolution_clock::now();

    int largeCoordSumCount = 0;

    const double startParam = nurbs.knotVector[0];
    const double endParam = nurbs.knotVector[nurbs.knotVector.size() - 1];

    std::vector<double> paramValues;

    // build a collection of parameter values which we will evaluate
    for (int i = 0; i < numSamples; i++) {
      const double t = startParam + i * (endParam - startParam) / (numSamples - 1);
      paramValues.push_back(t);
    }

    for (int i = 0; i < numEvals; i++) {
      // do a random evaluation
      const double t = paramValues[static_cast<int>(static_cast<double>(rand()) / RAND_MAX * numSamples)];
      const Point& result = nurbs.evaluate(t);

      // Prevent result from being optimized away
      (void)result;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    const double elapsedTime = duration.count();
    return elapsedTime;
  }

  const int NUM_SAMPLES = 50;
  const int NUM_EVALS = 20000;
  void doWork0() {
    const NURBS nurbs = degree2_nurbs_example1(0);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "No cache time\t" << elapsedTime << std::endl;
  }

  void doWork1() {
    // Plan: pass in a different configuration and expect different profiling results
    const NURBS nurbs = degree2_nurbs_example1(1);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "Map cache time \t" << elapsedTime << std::endl;
  }

  void doWork2() {
    const NURBS nurbs = degree2_nurbs_example1(2);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "Unordered map\t" << elapsedTime << std::endl;
  }

  void doWork3() {
    const NURBS nurbs = degree4_nurbs_example1(0);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "No cache time\t" << elapsedTime << std::endl;
  }

  void doWork4() {
    // Plan: pass in a different configuration and expect different profiling results
    const NURBS nurbs = degree4_nurbs_example1(1);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "Map cache time \t" << elapsedTime << std::endl;
  }

  void doWork5() {
    const NURBS nurbs = degree4_nurbs_example1(2);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "Unordered map\t" << elapsedTime << std::endl;
  }

  void doWork6() {
    const NURBS nurbs = degree6_nurbs_example1(0);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "No cache time\t" << elapsedTime << std::endl;
  }

  void doWork7() {
    // Plan: pass in a different configuration and expect different profiling results
    const NURBS nurbs = degree6_nurbs_example1(1);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "Map cache time \t" << elapsedTime << std::endl;
  }

  void doWork8() {
    const NURBS nurbs = degree6_nurbs_example1(2);
    const double elapsedTime = doConfiguredWork(nurbs, NUM_SAMPLES, NUM_EVALS);
    std::cout << __func__ << " " << "Unordered map\t" << elapsedTime << std::endl;
  }

  void nurbs_performance_example() {
    std::cout << "Degree 2 examples are quite fast to calculate basis values" << std::endl;
    doWork0();
    doWork1();
    doWork2();
    std::cout << "Degree 4 examples are slower to calculate basis values" << std::endl;
    doWork3();
    doWork4();
    doWork5();
    std::cout << "Degree 6 examples take more time to find basis values" << std::endl;
    doWork6();
    doWork7();
    doWork8();
  }

}