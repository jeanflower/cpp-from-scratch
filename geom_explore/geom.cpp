#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <tuple>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include "geom.hpp"

namespace geom_examples {
  // Define a 3D point structure
  Point::Point(double x, double y, double z) : x(x), y(y), z(z) {}

  double Point::X() const { return x; }
  double Point::Y() const { return y; }
  double Point::Z() const { return z; }

  Vector::Vector(double x, double y, double z) : x(x), y(y), z(z) {}

  double Vector::X() const { return x; }
  double Vector::Y() const { return y; }
  double Vector::Z() const { return z; }

  Circle::Circle(double r) : radius(r) {}

  Point Circle::evaluate(
    double t, 
    std::optional<Vector> first_derivative, 
    std::optional<Vector> second_derivative
  ) const {
    Point p{radius * cos(t), radius * sin(t)};

    if (first_derivative) {
      first_derivative -> x = -radius * sin(t);
      first_derivative -> y = radius * cos(t);
    }

    if (second_derivative) {
      second_derivative -> x = -radius * cos(t);
      second_derivative -> y = -radius * sin(t);
    }

    return p;
  }

  std::vector<geom_examples::PtCollection> viewables;

  void addGeometryToView(const std::vector<PtCollection>& ptColls) {
    for (auto x : ptColls) {
      viewables.push_back(x);
    }
  }

  void writeGeometryToJSON(
  ) {
    try {
      // Write this data to a json file for viewing to pick up
      std::ofstream debugFile("viewer/public/output/view_data.json");
      if (!debugFile) {
          std::cerr << "Error opening file for writing.\n";
          return;
      }

      std::vector<PtCollection> nonLineColls;
      std::vector<PtCollection> isLineColls;
      std::copy_if(
        viewables.begin(), 
        viewables.end(), 
        std::back_inserter(nonLineColls), 
        [](const PtCollection& elt) { return !elt.isLine; }
      );
      std::copy_if(
        viewables.begin(), 
        viewables.end(), 
        std::back_inserter(isLineColls), 
        [](const PtCollection& elt) { return elt.isLine; }
      );

      // we'll redirect cout to write to json file
      std::streambuf* coutBuf = std::cout.rdbuf(); // Save the original buffer

      std::cout.rdbuf(debugFile.rdbuf());  // Redirect std::cout to the file
      std::cout << "{";
      std::cout << "\"ptsGps\":[\n";
      for (size_t i = 0; i < nonLineColls.size(); ++i) {
        auto pts = nonLineColls[i].pts;
        auto col = nonLineColls[i].color;
        std::cout << "{";
        std::cout << "\"color\":" << col << ",\n";
        std::cout << "\"pts\":[\n";
        for (size_t j = 0; j < pts.size(); ++j) {
          const Point& pt = pts[j];
          std::cout 
            << "  { \"x\": " << pt.X() << ", \"y\": " << pt.Y() << ", \"z\": " << pt.Z() << "}";
          if (j != pts.size() - 1) {
            std::cout << ",";
          }
          std::cout << "\n"; // new point on a new line of the json file
        }
        std::cout << "]\n"; // end of array of points
        std::cout << "}\n";  // end of pts object
        if (i != nonLineColls.size() - 1) {
          std::cout << ",";
        }
      }
      std::cout << "],\n"; // end of array of pts objects
      std::cout << "\"linesObjs\":[\n";

      for (size_t polyLineIndex = 0; polyLineIndex < isLineColls.size(); ++polyLineIndex) {
        auto lineColl = isLineColls[polyLineIndex];
        auto col = lineColl.color;
        auto vxs = lineColl.pts;
        std::cout << "{";
        std::cout << "\"color\":" << col << ",\n";
        std::cout << "\"vxs\":";
        std::cout << "[\n";
        for (size_t segIndex = 0; segIndex < vxs.size(); ++segIndex) {
          const Point& pt = vxs[segIndex];
          std::cout 
            << "  { \"x\": " << pt.X() << ", \"y\": " << pt.Y() << ", \"z\": " << pt.Z() << "}";
          if (segIndex != vxs.size() - 1) {
            std::cout << ",";
          }
          std::cout << "\n"; // new point on a new line of the json file
        }
        std::cout << "]\n"; // end of vxs array
        std::cout << "}\n"; // end of polyLine object

        if (polyLineIndex != isLineColls.size() - 1) {
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
  double basisFunctionNoCache(const std::vector<double>& knot_vector, int i, int p, double t) {
    if (p == 0) { // zero degree basis function
      // Base case: the 0th degree basis function is 1 if t is within the interval [t_i, t_{i+1}]
      const double result = (t >= knot_vector[i] && t < knot_vector[i + 1]) ? 1.0 : 0.0;
      //indent(p);
      //std::cout << "N(" << i << ", " << p << ", " << t << ") = " << result << std::endl;;
      return result;
    } else { // higher degree basis function
      // Recursive case
      //indent(p);
      //std::cout << "N(" << i << ", " << p << ", " << t << ") needs lower degree work ..." << std::endl;

      double denom1 = knot_vector[i + p] - knot_vector[i];
      double denom2 = knot_vector[i + p + 1] - knot_vector[i + 1];

      double num1 = (t - knot_vector[i]);
      double num2 = (knot_vector[i + p + 1] - t);

      double basis1 = basisFunctionNoCache(knot_vector, i, p - 1, t);
      double basis2 = basisFunctionNoCache(knot_vector, i + 1, p - 1, t);

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

  enum CacheStrategy {
    NO_CACHE = 0,
    MAP_CACHE_SIMPLE = 1,
    MAP_CACHE_EMPLACE = 2,
    UNORDERED_MAP_CACHE = 3
  };

  class BasisEvaluator {
    private:
      std::vector<double> knot_vector;

    public:
      int cache_strategy = NO_CACHE;
      // strategy 0 uses no cache
      std::map<std::tuple<int, int, double>, double> basisCacheMap; 
      std::unordered_map<std::tuple<int, int, double>, double, KeyHash> basisCacheUnorderedMap;

      BasisEvaluator(
        const std::vector<double>& knot_vector,
        int cache_strategy = NO_CACHE
      ): knot_vector(knot_vector), cache_strategy(cache_strategy) {}

      double basisFunction(int i, int p, double t) {
        if (cache_strategy == NO_CACHE) {
          return basisFunctionNoCache(knot_vector, i, p, t);
        } else if(cache_strategy == MAP_CACHE_SIMPLE) {
          // Use a std::map as a cache to store the basis function values
          auto key = std::make_tuple(i, p, t);
          if (basisCacheMap.find(key) == basisCacheMap.end()) {
            // std::cout << "Cache miss for basis function (" << i << ", " << p << ", " << t << ")" << std::endl;
            basisCacheMap[key] = basisFunctionNoCache(knot_vector, i, p, t);
          } else {
            // std::cout << "Cache hit for basis function (" << i << ", " << p << ", " << t << ")" << std::endl;
          }
          return basisCacheMap[key];    
        } else if(cache_strategy == MAP_CACHE_EMPLACE) {
          // Use a std::map as a cache to store the basis function values
          auto key = std::make_tuple(i, p, t);

          // Use emplace to avoid redundant lookups
          auto [iter, inserted] = basisCacheMap.emplace(key, 0.0); // Insert a default value (0.0) if key is missing
          if (inserted) {
              iter->second = basisFunctionNoCache(knot_vector, i, p, t);
          }
          return iter->second;
        } else if (cache_strategy == UNORDERED_MAP_CACHE) {
            // Use an unordered_map-based cache for faster lookups
            auto key = std::make_tuple(i, p, t);
            auto [iter, inserted] = basisCacheUnorderedMap.emplace(key, 0.0); // Insert a default value (0.0) if key is missing
            if (inserted) {
                iter->second = basisFunctionNoCache(knot_vector, i, p, t);
            }
            return iter->second;
          } else {
          return basisFunctionNoCache(knot_vector, i, p, t);
        }
      }
  };

  // Define the NURBS curve class
  class NURBS {
    public:
      int degree;
      std::vector<Point> control_points;
      std::vector<double> weights;
      std::vector<double> knot_vector;

    private:
      bool useCache = false;
      mutable BasisEvaluator basisEvaluator;

    public:
      NURBS(
        int degree, 
        const std::vector<Point>& control_points, 
        const std::vector<double>& weights, 
        const std::vector<double>& knot_vector,
        int cache_strategy
      ): degree(degree), control_points(control_points), 
          weights(weights), knot_vector(knot_vector), 
          basisEvaluator(knot_vector, cache_strategy) {
        if (degree + control_points.size() + 1 != knot_vector.size()) {
          std::cerr << "Invalid NURBS curve definition: degree + control points + 1 != knot vector size" << std::endl;
        }
        if (weights.size() != control_points.size()) {
          std::cerr << "Invalid NURBS curve definition: weights size != control points size" << std::endl;
        }
      }

      double basisFunction(int i, int p, double t) const {
        return basisEvaluator.basisFunction(i, p, t);
      }

      // Evaluate the NURBS curve at parameter t
      const Point evaluate(double t) const {
        int n = control_points.size();
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
          numeratorX += N * weight * control_points[i].x;
          numeratorY += N * weight * control_points[i].y;
          numeratorZ += N * weight * control_points[i].z;
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

  NURBS degree2NurbsExample1(
    int cache_strategy
  ) {
    // Define degree, control points, weights, and knot vector
    int degree = 2;
    std::vector<Point> control_points = {
      Point(0, 0, 0), 
      Point(1, 1, 0), 
      Point(2, 0, 0)
    };
    std::vector<double> weights = {
      1.0, 
      1.0, 
      1.0
    };
    std::vector<double> knot_vector = {  // 3 + 2 + 1 = 6 knots
      0.0,
      0.0,
      0.0,
      // 0.5, too many knots! Fails reading control points out-of-range
      1.0,
      1.0,
      1.0
    };
    BasisEvaluator basisEvaluator(knot_vector, cache_strategy);
    // Create the NURBS curve
    NURBS nurbs(degree, control_points, weights, knot_vector, cache_strategy);
    return nurbs;
  }

  NURBS degree3NurbsExample1(
    int cache_strategy
  ) {
    // Define degree, control points, weights, and knot vector
    int degree = 3;
    std::vector<Point> control_points = {
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
    std::vector<double> knot_vector = {   // 4 + 3 + 1 = 8 knots
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
    NURBS nurbs(degree, control_points, weights, knot_vector, cache_strategy);
    return nurbs;
  }

  NURBS degree4NurbsExample1(
    int cache_strategy
  ) {
    // Define degree, control points, weights, and knot vector
    int degree = 4;
    std::vector<Point> control_points = {
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
    std::vector<double> knot_vector = { // 5 + 4 + 1 = 10 knots
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
    NURBS nurbs(degree, control_points, weights, knot_vector, cache_strategy);
    return nurbs;
  }

  NURBS degree6NurbsExample1(
    int cache_strategy
  ) {
    // Define degree, control points, weights, and knot vector
    int degree = 6;
    std::vector<Point> control_points = {
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
    std::vector<double> knot_vector = { // 7 + 6 + 1 = 14 knots
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
    NURBS nurbs(degree, control_points, weights, knot_vector, cache_strategy);
    return nurbs;
  }

  void nurbsExample() {

    /*
    // check for read out of bounds
    std::vector<double> weights = {
      1.0, 
      1.0, 
      1.0
    };
    std::cout << "weights[3] = " << weights[3] << std::endl;
    */
   const NURBS nurbs = degree2NurbsExample1(NO_CACHE);

    /*
    // Basis functions should sum to 1
    const double t = 0.499;
    double sum = 0.0;
    for (int i = 0; i < nurbs.control_points.size(); ++i) {
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
    const double start = nurbs.knot_vector[0];
    const double end = nurbs.knot_vector[nurbs.knot_vector.size() - 1];
    for (int i = 0; i < NUM_SAMPLES; i++) {
      const double t = start + i * (end - start) / (NUM_SAMPLES - 1);
      const Point& result = nurbs.evaluate(t);
      ptsToDisplay.push_back(Point(result.x, result.y, result.z));
      // std::cout << t << ", " 
      // << result.x << ", " << result.y << ", " << result.z << std::endl;
    }
    std::vector<PtCollection> ptsColls;
    PtCollection ptsObj ={
      .color = GREEN,
      .pts = ptsToDisplay,
      .isLine = false
    };
    ptsColls.push_back(ptsObj);
    addGeometryToView(ptsColls);
  }

  double doConfiguredWork(
    const NURBS& nurbs,
    const int numParamValsPossible, 
    const int numEvals
  ) {
    auto start = std::chrono::high_resolution_clock::now();

    int largeCoordSumCount = 0;

    const double startParam = nurbs.knot_vector[0];
    const double endParam = nurbs.knot_vector[nurbs.knot_vector.size() - 1];

    std::vector<double> paramValues;

    // build a collection of parameter values which we will evaluate
    for (int i = 0; i < numParamValsPossible; i++) {
      const double t = startParam + i * (endParam - startParam) / (numParamValsPossible - 1);
      paramValues.push_back(t);
    }

    for (int i = 0; i < numEvals; i++) {
      // do a random evaluation
      const double t = paramValues[static_cast<int>(static_cast<double>(rand()) / RAND_MAX * numParamValsPossible)];
      const Point& result = nurbs.evaluate(t);

      // Prevent result from being optimized away
      (void)result;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    const double elapsedTime = duration.count();
    return elapsedTime;
  }

  std::string noCacheMessage(double t) {
    return "  No cache time            " + std::to_string(t) + "\n";
  }
  std::string mapSimpleCacheMessage(double t) {
    return "  Simple map cache time    " + std::to_string(t) + "\n";
  }
  std::string mapEmplaceCacheMessage(double t) {
    return "  Emplace map cache time   " + std::to_string(t) + "\n";
  }
  std::string unorderedMapCacheMessage(double t) {
    return "  Unordered map cache time " + std::to_string(t) + "\n";
  }
  std::string doWork(int d, int cache_strategy) {
    const NURBS nurbs = 
      d == 2 ? degree2NurbsExample1(cache_strategy) :
      d == 4 ? degree4NurbsExample1(cache_strategy) :
      degree6NurbsExample1(cache_strategy);

    // a big number of evaluations
    //  - uses the cache more
    //  - reduces the effect of noise in the performance measures
    const int numEvaluationsForPerf = 20000;

    // a limited number of possible parameter values
    //  - uses the cache more
    //  - controls the size of the cache
    //  - makes the test less realistic for some solving scenarios
    const int numParamValsPossible = 1000;
    const double elapsedTime = doConfiguredWork(nurbs, numParamValsPossible, numEvaluationsForPerf);

    std::string result = "";
    if (cache_strategy == NO_CACHE) {
      result += noCacheMessage(elapsedTime);
    } else if (cache_strategy == MAP_CACHE_SIMPLE) {
      result += mapSimpleCacheMessage(elapsedTime);
    } else if (cache_strategy == MAP_CACHE_EMPLACE) {
      result += mapEmplaceCacheMessage(elapsedTime);
    } else if (cache_strategy == UNORDERED_MAP_CACHE) {
      result +=  unorderedMapCacheMessage(elapsedTime);
    }
    std::cout << result;
    return result;
  }

  std::string nurbsPerformanceExample() {

    std::string result = "";
    std::string heading = "Degree 2 examples are quite fast to calculate basis values\n";
    std::cout << heading;
    result += heading;
    result += " " + doWork(2, NO_CACHE);
    result += " " + doWork(2, MAP_CACHE_SIMPLE);
    result += " " + doWork(2, MAP_CACHE_EMPLACE);
    result += " " + doWork(2, UNORDERED_MAP_CACHE);
    heading = "Degree 4 examples are slower to calculate basis values\n";
    std::cout << heading;
    result += heading;
    result += " " + doWork(4, NO_CACHE);
    result += " " + doWork(4, MAP_CACHE_SIMPLE);
    result += " " + doWork(4, MAP_CACHE_EMPLACE);
    result += " " + doWork(4, UNORDERED_MAP_CACHE);
    heading = "Degree 6 examples take more time to find basis values\n";
    std::cout << heading;
    result += heading;
    result += " " + doWork(6, NO_CACHE);
    result += " " + doWork(6, MAP_CACHE_SIMPLE);
    result += " " + doWork(6, MAP_CACHE_EMPLACE);
    result += " " + doWork(6, UNORDERED_MAP_CACHE);

    return result;
  }

  void circleExample() {
    Circle c(5.0);
    double t = 1.0;

    Vector first, second;
    Point p = c.evaluate(t, first, second);

    std::cout << "Point: (" << p.x << ", " << p.y << ")\n";
    std::cout << "First Derivative: (" << first.x << ", " << first.y << ")\n";
    std::cout << "Second Derivative: (" << second.x << ", " << second.y << ")\n";
  }
}