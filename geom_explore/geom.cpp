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
#include "geom_display.hpp"
#include <functional>
#include <utility>  // for std::pair
#include <string>
#include <vector>
#include <mutex>
#include <thread>

namespace geom_examples {
  // Define a 3D point structure
  Point::Point(double x, double y, double z) : x(x), y(y), z(z) {}

  double Point::X() const { return x; }
  double Point::Y() const { return y; }
  double Point::Z() const { return z; }

  // Define a 3D vector structure - same as a Point in terms of data but
  // in a context, means something different
  Vector::Vector(double x, double y, double z) : x(x), y(y), z(z) {}

  double Vector::X() const { return x; }
  double Vector::Y() const { return y; }
  double Vector::Z() const { return z; }

  // Define a Line in 3d space from a start point
  // parameterised along a given direction
  // the length of the direction vector provides the parameterisation
  Line::Line(
    Point start,
    Vector direction
  ) : start(start), direction(direction) {}

  // A Line can be evaluated and provide its derivatives
  Point Line::evaluate(
    double t, 
    std::optional<std::reference_wrapper<Vector>> first_derivative, 
    std::optional<std::reference_wrapper<Vector>> second_derivative
  ) const {
    Point p{
      start.X() + t * direction.X(), 
      start.Y() + t * direction.Y(), 
      start.Z() + t * direction.Z() 
    };

    if (first_derivative.has_value()) {
      first_derivative->get().x = direction.X();
      first_derivative->get().y = direction.Y();
      first_derivative->get().z = direction.Z();
    }

    if (second_derivative.has_value()) {
      second_derivative->get().x = 0;
      second_derivative->get().y = 0;
      second_derivative->get().z = 0;
    }

    // deliberately break the evaluator to test our tests!!
    //if (t == 1000.0) { // terrible way to pick out a t
    //  p.x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && first_derivative.has_value()) { // terrible way to pick out a t
    //  first_derivative->get().x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && second_derivative.has_value()) { // terrible way to pick out a t
    //  second_derivative->get().x = 88; // nonsense value to return
    //}
    return p;
  }

  // A Circle curve centered at the given center, with given radius
  // it always falls in the x/y plane and is parameterised in the same way
  Circle::Circle(
    Point center,
    double r
  ) : center(center), radius(r) {}

  // A Circle can be evaluated and provide its derivatives
  Point Circle::evaluate(
    double t, 
    std::optional<std::reference_wrapper<Vector>> first_derivative, 
    std::optional<std::reference_wrapper<Vector>> second_derivative
  ) const {
    Point p{
      center.X() + radius * cos(t), 
      center.Y() + radius * sin(t),
      center.Z()
    };

    if (first_derivative.has_value()) {
      first_derivative->get().x = -radius * sin(t);
      first_derivative->get().y = radius * cos(t);
    }

    if (second_derivative.has_value()) {
      second_derivative->get().x = -radius * cos(t);
      second_derivative->get().y = -radius * sin(t);
    }

    // deliberately break the evaluator to test our tests!!
    //if (t == 1000.0) { // terrible way to pick out a t
    //  p.x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && first_derivative.has_value()) { // terrible way to pick out a t
    //  first_derivative->get().x = 88; // nonsense value to return
    //}
    //if (t == 1000.0 && second_derivative.has_value()) { // terrible way to pick out a t
    //  second_derivative->get().x = 88; // nonsense value to return
    //}
    return p;
  }
  
  // Build a circle object and evaluate it
  void circleExample() {
    Circle c(Point(0.0, 0.0, 0.0), 5.0);
    double t = 1.0;

    Vector first, second;
    Point p = c.evaluate(t, first, second);

    std::cout << "Point: (" << p.x << ", " << p.y << ")\n";
    std::cout << "First Derivative: (" << first.x << ", " << first.y << ")\n";
    std::cout << "Second Derivative: (" << second.x << ", " << second.y << ")\n";
  }

  // Compare the error between the geometry evaluated at t+dt
  // and the prediction based on the geometry at t and the derivative
  // giving an approximation of the geometry at t+dt
  // This function is used to validate the evaluated derivatives
  double calculate_sq_error(
    std::function<std::pair<Vector, Vector>(double)>& geom_evaluator,
    double t,
    const Vector p,
    const Vector& first_deriv,
    double t_step
  ){
    std::pair<Vector, Vector> evals = geom_evaluator(t + t_step); // above t
    Vector p_nearby = evals.first; // reads like position but could actually be a higher deriv
    // std::cout << "p_nearby = " << p_nearby.X() << ", " << p_nearby.Y() << ", " << p_nearby.Z() << "\n";

    Vector estimate( // TODO implement subtraction of Point
      (p_nearby.X() - p.X()) / t_step,
      (p_nearby.Y() - p.Y()) / t_step,
      (p_nearby.Z() - p.Z()) / t_step
    );

    const double sq_error =  // TODO implement subtraction, sq_length of Vector
      (estimate.X() - first_deriv.X()) * (estimate.X() - first_deriv.X()) +
      (estimate.Y() - first_deriv.Y()) * (estimate.Y() - first_deriv.Y()) +
      (estimate.Z() - first_deriv.Z()) * (estimate.Z() - first_deriv.Z());

    std::cout << "estimate with t_step " << t_step << " = "
      << estimate.X() << ", " << estimate.Y() << ", " << estimate.Z()
      << " with sq error " << sq_error
      << "\n";
    return sq_error;
  }

  // Can be used for testing derivatives at any order
  // being given evaluators for order n and n+1
  bool test_curve_derivs_generic(
    std::function<std::pair<Vector, Vector>(double)>& geom_evaluator,
    double t
  ) {
    std::pair<Vector, Vector> evals = geom_evaluator(t);
    Vector p = evals.first; // reads like position but could actually be a higher deriv
    Vector first_deriv = evals.second;

    //std::cout << "p = " << p.X() << ", " << p.Y() << ", " << p.Z() << "\n";
    //std::cout << "first_deriv = " << first_deriv.X() << ", " << first_deriv.Y() << ", " << first_deriv.Z() << "\n";

    const double sq_epsilon = 1e-9;

    // convergence means: for all epsilon, there exists a delta such that...
    // for all values of t within delta
    // the error between estimated deriv and reported deriv is less than epsilon
    //
    // we're just picking one epsilon, guessing a delta,
    // and performing a sample of evaluations within delta
    // this isn't quite right in many ways...
    // it assumes a monotonic convergence from our first estimate within tolerance

    // ensure that, if we get close enough, the estimate is close to the reported value
    // to this tolerance
    const double minimal_delta = 1e-8; // if epsilon gets this small, we are close
    const double initial_t_step = 0.001; // start relatively close for a good chance of a good estimate
    const double initial_large_sq_error = std::numeric_limits<double>::max();
    const double t_step_scaling = 0.5; // reduce t_step each time

    // use a parameter-step of t_step to estimate the derivative at t
    double t_step = initial_t_step; 

    double sq_error = initial_large_sq_error;
    bool found_good_estimate = false;
    while (sq_error > sq_epsilon && t_step > minimal_delta) {
      t_step *= t_step_scaling;
      sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      found_good_estimate = sq_error < sq_epsilon;
    }

    if (!found_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, no estimate from above\n";
      return false;
    }

    //std::cout << "found good estimate with delta = " << t_step << "\n";
    bool stayed_good_estimate = true;
    while (stayed_good_estimate && t_step > minimal_delta) {
      t_step *= t_step_scaling;
      const double sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      stayed_good_estimate = sq_error < sq_epsilon;
    }
    if (!stayed_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, lost estimate from above\n";
      return false;
    }
    //std::cout << "kept good estimates to delta = " << t_step << "\n";

    t_step = -initial_t_step; 

    sq_error = initial_large_sq_error;
    found_good_estimate = false;
    while (sq_error > sq_epsilon && -t_step > minimal_delta) {
      t_step *= t_step_scaling;
      sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      found_good_estimate = sq_error < sq_epsilon;
    }

    if (!found_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, no estimate from below\n";
      return false;
    }
    //std::cout << "found good estimate with t_step = " << t_step << "\n";

    stayed_good_estimate = true;
    while (stayed_good_estimate && -t_step > minimal_delta) {
      t_step *= t_step_scaling;
      const double sq_error = calculate_sq_error(
        geom_evaluator, t, p, first_deriv, t_step
      );
      stayed_good_estimate = sq_error < sq_epsilon;
    }
    if (!stayed_good_estimate) {
      std::cout << "Failed to converge towards first_deriv, lost estimate from below\n";
      return false;
    }
    //std::cout << "kept good estimates to delta = " << t_step << "\n";

    return true;
  }

  bool test_curve_derivs_at(
    const Curve& c,
    const double t
  ) {
    bool result = true;
    std::cout << "Test curve derivs at parameter value " << t << "\n";

    std::function<std::pair<Vector, Vector>(double)> position_and_first_deriv_evaluator =
      [&c](double t) -> std::pair<Vector, Vector> {
      Vector firstDeriv;
      Point p = c.evaluate(t, firstDeriv);
      const Vector v(p.X(), p.Y(), p.Z()); // TODO avoid ever copying Point to Vector
      return std::make_pair(
        v, 
        firstDeriv
      );
    };
    result = test_curve_derivs_generic(
      position_and_first_deriv_evaluator,
      t
    );

    if (result) {
      std::function<std::pair<Vector, Vector>(double)> first_and_second_deriv_evaluator =
        [&c](double t) -> std::pair<Vector, Vector> {
        Vector firstDeriv, secondDeriv;
        Point p = c.evaluate(t, firstDeriv, secondDeriv);
        return std::make_pair(
          firstDeriv,
          secondDeriv
        );
      };
      result = test_curve_derivs_generic(
        first_and_second_deriv_evaluator,
        t
      );
    }
    return result;
  }

  // We'd expect this to pass unless we 'nobble' the Circle evaluators
  // to intentionally return bad values for specific parameter values
  bool test_curve_derivs(){

    // select a curve
    const Circle c(Point(0.0, 0.0, 0.0), 12.0);

    // select values of t to analyse for derivative accuracy    
    bool result = true;
    if (result) {
      result = test_curve_derivs_at(c, 0.2);
    }
    if (result) {
      result = test_curve_derivs_at(c, 0);
    }
    if (result) {
      result = test_curve_derivs_at(c, 1000.0);
    }

    return result;
  }

  // We might seek the closest point to the origin...
  // or solve for where a curve passes through the origin
  // this creates a function of one variable returning one number
  double DistanceSqFromOrigin::operator() (const double t) {
    Point p = c.evaluate(t);
    return 
      p.X() * p.X() +
      p.Y() * p.Y() +
      p.Z() * p.Z();
  };
  double DistanceSqFromOrigin::df(const double t) {
    Vector deriv;
    Point p = c.evaluate(t, deriv);
    return
      2 * p.X() * deriv.X() +
      2 * p.Y() * deriv.Y() +
      2 * p.Z() * deriv.Z();
  };

  // Newton Raphson seeks a root of a function of one variable
  template <class T>
  double rtnewt(
    T &funcd, 
    const double guess,
    const double lower_bound_param, 
    const double upper_bound_param, 
    const double accuracy_tolerance
  ) {
    const int MAX_ITERATIONS = 200;
    // initial guess
    double rtn = guess;
    for (int j = 0; j < MAX_ITERATIONS; j++) {
      double f = funcd(rtn);
      double df = funcd.df(rtn);
      if (df == 0.0) 
        throw std::runtime_error("df == 0.0 in rtnewt");
      double dx = f / df;
      rtn -= dx;
      // std::cout << "j is " << j << " and rtn is " << rtn << "\n";
      if ((lower_bound_param - rtn) * (rtn - upper_bound_param) < 0.0)
        throw std::runtime_error("Jumped out of range in rtnewt after " + std::to_string(j) + " iterations");
      if (abs(dx) < accuracy_tolerance) 
        return rtn; 
    }
    throw("Maximum number of iterations exceeded in rtnewt");
  }

  // This function builds some curves and uses Newton Raphson to seek
  // a point where those curves pass through the origin.
  // Examples may or may not converge.
  void test_distance_sq_fn() {

    std::function<void(DistanceSqFromOrigin&, double)> print_vals =
      [](DistanceSqFromOrigin& f, double t) -> void {
        const double val = f(t);
        const double dval = f.df(t);
    
        std::cout << "distance from O at " << t << " is " << val << " and df is " << dval << "\n";
      };

    Circle c1(Point(0.0, 0.0, 0.0), 2);
    DistanceSqFromOrigin f1(c1);
  
    //print_vals(f1, 0);
    //print_vals(f1, 3);

    Circle c2(Point(1.0, 2.0, 0.0), sqrt(5.0));
    DistanceSqFromOrigin f2(c2);
  
    //print_vals(f2, 0);
    //print_vals(f2, 3);

    Circle c3(Point(1.0, 2.0, 3.0), sqrt(14.0));
    DistanceSqFromOrigin f3(c3);
  
    //print_vals(f2, 0);
    //print_vals(f2, 3);

    addCurveToView(Line(Point(0,0,0), Vector(1,0,0)), 0.0, 1.0, RED, 1);
    addCurveToView(Line(Point(0,0,0), Vector(0,1,0)), 0.0, 1.0, GREEN, 1);
    addCurveToView(Line(Point(0,0,0), Vector(0,0,1)), 0.0, 1.0, BLUE, 1);

    addCurveToView(c1, 0.0, 2*M_PI, GREEN, 2);
    addCurveToView(c2, 0.0, 2*M_PI, RED, 2);
    addCurveToView(c3, 0.0, 2*M_PI, BLUE, 2);

    // e.g. we might highlight a 'closest point'
    addPointToView(Point(2, 0, 0), YELLOW, 7);

    // f1 doesn't have a zero - don't expect NR to converge
    try {
      double newtonResult = rtnewt(f1, M_PI, -100.0 * M_PI, 100.0 * M_PI, 1e-6);
      std::cout << "Unepected newtonResult = " << newtonResult << "\n";
    } catch (const std::runtime_error& e) {
      std::cerr << "As expeced for f1, get an error: " << e.what() << std::endl;
    } catch (const char* msg) {
      std::cerr << "As expeced for f1, get an error: " << msg << std::endl;
    } catch (...) {
      std::cerr << "Unknown exception caught!" << std::endl;
    }

    // f2 does have a zero - expect NR to converge
    for (int j = 0; j < 10; j++) {
      try {
        double guess = j * 0.1 * 2 * M_PI;
        std::cout << "newtonResult on f2 from guess " << guess << " ";
        double newtonResult = rtnewt(f2, guess, -100.0 * M_PI, 100.0 * M_PI, 1e-6);
        std::cout << newtonResult << "\n";
      } catch (const std::runtime_error& e) {
        std::cerr << "Unexpected for f2, get an error: " << e.what() << std::endl;
        } catch (const char* msg) {
        std::cerr << "Unexpected for f2, error: " << msg << std::endl;
      } catch (...) {
        std::cerr << "Unknown exception caught!" << std::endl;
      }
    }

    geom_examples::writeGeometryToJSON();
  }

  // A function of two variables uses (x,y)
  // and a function of a complex variable uses x+iy
  template <class T>
  struct Coords2D {
    T x, y;

    // Constructor
    Coords2D(T x = 0, T y = 0);

    // Getter functions
    T X() const;
    T Y() const;
  };

  template <class T>
  Coords2D<T>::Coords2D(T x, T y) : x(x), y(y) {}

  template <class T>
  T Coords2D<T>::X() const { return x; }
  template <class T>
  T Coords2D<T>::Y() const { return y; }

  template <class T>
  Coords2D<T> operator+(const Coords2D<T>& a, const Coords2D<T>& b) {
    return Coords2D(a.x + b.x, a.y + b.y);
  }

  template <class T>
  Coords2D<T> operator-(const Coords2D<T>& a, const Coords2D<T>& b) {
    return Coords2D(a.x - b.x, a.y - b.y);
  }

  // A certain kind of cubic function of two real values
  // inspired by the complex function z^3 - 1
  template <class T>
  class CubicFunction {
    public:
      explicit CubicFunction(
        std::vector<T> coeffs,
        std::vector<Coords2D<T>> roots
      ) {
        this->roots = roots;
        // ax^3 + bx^2y + cxy^2 + dy^3 + e
        // fx^3 + gx^2y + hxy^2 + iy^3 + j
        a = coeffs[0]; 
        b = coeffs[1]; 
        c = coeffs[2];
        d = coeffs[3];
        e = coeffs[4];
        f = coeffs[5];
        g = coeffs[6];
        h = coeffs[7];
        i = coeffs[8];
        j = coeffs[9];
      }

      std::vector<Coords2D<T>> roots;

      private: 
        // coefficients of the cubic polynomial
        T a = 0.0, b = 0.0, c = 0.0, d = 0.0;
        T e = 0.0;

        T f = 0.0, g = 0.0, h = 0.0, i = 0.0;
        T j = 0.0;

      public: 

        // expand a cubic polynomial
        Coords2D<T> operator() (const Coords2D<T>& z) {
          return Coords2D<T>(          
              a * z.X() * z.X() * z.X() // a = 1    1 * z.X() * z.X() * z.X()
            + b * z.X() * z.X() * z.Y() // b = 0    0
            + c * z.X() * z.Y() * z.Y() // c = -3  -3 * z.X() * z.Y() * z.Y()
            + d * z.Y() * z.Y() * z.Y() // d = 0    0
            + e,                        // e = -1  -1
              f * z.X() * z.X() * z.X() // f = 0    0
            + g * z.X() * z.X() * z.Y() // g = 3    3 * z.X() * z.X() * z.Y()
            + h * z.X() * z.Y() * z.Y() // h = 0    0
            + i * z.Y() * z.Y() * z.Y() // i = -1   -1 * z.Y() * z.Y() * z.Y()
            + j                         // j = 0    0
          );
        }
        // derivative of a cubic polynomial
        Coords2D<T> dfx(const Coords2D<T>& z) {
          return Coords2D<T>(
              3 * a * z.X() * z.X() // a = 1    1 * 3 * z.X() * z.X()
            + 2 * b * z.X() * z.Y() // b = 0    0
            +     c * z.Y() * z.Y() // c = -3   -3 * z.Y() * z.Y()
            ,
              3 * f * z.X() * z.X() // f = 0    0
            + 2 * g * z.X() * z.Y() // g = 3    2 * 3 * z.X() * z.Y()
            +     h * z.Y() * z.Y() // h = 0    0
          );
        }
        // derivative of a cubic polynomial
        Coords2D<T> dfy(const Coords2D<T>& z) {
          return Coords2D<T>(

                  b * z.X() * z.X() // b = 0    0
            + 2 * c * z.X() * z.Y() // c = -3   -6 * z.X() * z.Y()
            + 3 * d * z.Y() * z.Y() // d = 0    0
            ,
                  g * z.X() * z.X() // g = 3    3 * z.X() * z.X()
            + 2 * h * z.X() * z.Y() // h = 0    0
            + 3 * i * z.Y() * z.Y() // i = -1   -1 * 3 * z.Y() * z.Y()
          );
        }
  };

  // Seek a root of a function of two variables
  // Templatise so that we can explore floats, doubles, long doubles
  template <class T, class F>
  Coords2D<T> newtonRaphson2Dinput(
    F &funcd, 
    Coords2D<T>& guess,
    const T lower_bound_x, 
    const T upper_bound_x, 
    const T lower_bound_y, 
    const T upper_bound_y, 
    const T accuracy_tolerance,
    const int MAX_ITERATIONS
  ) {
    const double accuracy_tolerance_sq = accuracy_tolerance * accuracy_tolerance;
    Coords2D<T> rtn = guess;
    for (int j = 0; j < MAX_ITERATIONS; j++) {
      Coords2D<T> f = funcd(rtn);
      Coords2D<T> dfx = funcd.dfx(rtn);
      Coords2D<T> dfy = funcd.dfy(rtn);

      double det = dfx.X() * dfy.Y() - dfx.Y() * dfy.X();
      if(det == 0) {
          throw std::runtime_error("Zero derivative determinant in complex_newt");
      }
      
      double step_x = ( f.X() * dfy.Y() + f.Y() * dfx.Y() ) / det;
      double step_y = ( f.Y() * dfx.X() + f.X() * dfy.X() ) / det;
      Coords2D<T> step(step_x, step_y);

      rtn = rtn - step;

      // std::cout << "j is " << j << " and rtn is " << rtn.X() << ", " << rtn.Y() << "\n";
      if ((lower_bound_x - rtn.X()) * (rtn.X() - upper_bound_x) < 0.0 ||
        (lower_bound_y - rtn.Y()) * (rtn.Y() - upper_bound_y) < 0.0)
        throw("Jumped out of range in newtonRaphson2Dinput");
      if (step.X() * step.X() + step.Y() * step.Y() < accuracy_tolerance_sq) {
        // std::cout << "Converged to " << rtn.X() << ", " << rtn.Y() << " after " << j << " iterations\n";
        return rtn; 
      }
    }
    throw("Maximum number of iterations exceeded in newtonRaphson2Dinput");
  }

  // Assess whether the 2D Newton Raphson converges for a given function f
  // starting from a given guess.
  // Return a number which says which of the knownSolutions convergence gave
  // or one more if we didn't converge to a known solution.
  // e..g. provide 3 known roots, and this can return 0, 1 or 2 for convergence to a known root
  // and 4 for convergence somewhere else or failure to converge.
  template <class T>
  int NR2DConverges(
    CubicFunction<T>& f,
    Coords2D<T>& guess,
    T accuracy_tolerance,
    const std::vector<Coords2D<T>>& knownSolutions
  ) {
    try {
      Coords2D<T> start = guess;

      //std::cout << "newtonResult on f from guess " << start.X() << " + i * " << start.Y() << "\n";
      Coords2D<T> newtonResult = newtonRaphson2Dinput<T>(
        f, 
        start,
        -100.0, 100.0, 
        -100.0, 100.0, 
        accuracy_tolerance,
        100
      );
      //std::cout << newtonResult.X() << " + i * " << newtonResult.Y() << "\n";

      for (int i = 0; i < knownSolutions.size(); i++) {
        const Coords2D<T>& knownSoln = knownSolutions[i];
        if (abs(newtonResult.X() - knownSoln.X()) < 0.1 
         && abs(newtonResult.Y() - knownSoln.Y()) < 1e-6) {
          // std::cout << "Converged to " << newtonResult.X() << " + i * " << newtonResult.Y() << " on root " << i << "\n";  
          // we converged on the root of interest
          return i;
        }
      }
      return knownSolutions.size();

    } catch (const char* msg) {
      //std::cerr << "Error: " << msg << std::endl;
    } catch (...) {
      //std::cerr << "Unknown exception caught!" << std::endl;
    }
    return false;
  }

  // Assess whether the 2D Newton Raphson converges for a given function f
  // starting from a range of guesses in the given low/high ranges.

  // Return a collection of Coords2D for guesses which converge to each of known solutions.
  // e..g. provide 3 known roots, and this can return 
  // - three collections of Coords2D which for convergence to a known root
  // - and 4th collection of Coords2D which converge somewhere else or fail to converge.
  template<class T>
  void assessConvergence(
    CubicFunction<T> f,
    std::vector<std::vector<Coords2D<T>>>& ptsData,
    int NUM_I, 
    int NUM_J, 
    T LOW_X, 
    T HIGH_X, 
    T LOW_Y, 
    T HIGH_Y,
    T accuracy_tolerance
  ) {
    std::mutex mtx; // protects access to ptsConverged and ptsNotConverged

    // Lambda that processes a block of rows (i values)
    auto process_range = [&](int start_i, int end_i) {
      std::vector<std::vector<Coords2D<T>>> localPtsData(4);

      for (int i = start_i; i < end_i; i++) {
        for (int j = 0; j < NUM_J; j++) {
          Coords2D<T> start(
            LOW_X + (HIGH_X - LOW_X) / NUM_I * i,
            LOW_Y + (HIGH_Y - LOW_Y) / NUM_J * j
          );
          int converged = NR2DConverges<T>(
            f, 
            start,
            accuracy_tolerance,
            f.roots
          );
          localPtsData[converged].push_back(start);
        }
      }
      // Merge local results into the global vectors with locking
      std::lock_guard<std::mutex> lock(mtx);
      for (int i = 0; i < 4; i++) {
        ptsData[i].insert(ptsData[i].end(), localPtsData[i].begin(), localPtsData[i].end());
      }
    };

    // Determine the number of threads available.
    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) { numThreads = 2; } // Fallback in case hardware_concurrency returns 0

    std::vector<std::thread> threads;
    int range = NUM_I / numThreads;
    for (unsigned int t = 0; t < numThreads; t++) {
      int start_i = t * range;
      int end_i = (t == numThreads - 1) ? NUM_I : start_i + range;
      threads.emplace_back(process_range, start_i, end_i);
    }
    // Join all threads.
    for (auto &th : threads) {
      th.join();
    }    
  }

  // Given collections of Points, assign different colors and add them to the viewables
  void addPtsToView(
    std::vector<std::vector<Point>>& ptsData,
    int displaySize
  ){
    std::vector<PtCollection> ptsColls;

    std::vector<uint32_t> colors = {
      geom_examples::RED, 
      geom_examples::GREEN, 
      geom_examples::BLUE, 
      geom_examples::YELLOW
    };
    for (int i = 0; i < 4; i++) {
      PtCollection ptsObj = {
        .displaySize = displaySize,
        .color = colors[i],
        .pts = ptsData[i],
        .isLine = false
      };
      ptsColls.push_back(ptsObj);
    }

    // std::cout << "adding " << ptsObj.pts.size() << " to view\n";
    addGeometryToView(ptsColls);    
  }

  // An example of a function of two variables
  template<class T>
  CubicFunction<T> zCubedMinus1 = CubicFunction<T>(
    // "z^3-1",
    { 
      1, 0, -3, 0, -1, 
      0, 3, 0, -1, 0
    },
    std::vector<Coords2D<T>>{
      Coords2D<T>(1, 0),
      Coords2D<T>(-0.5, sqrt(3) / 2),
      Coords2D<T>(-0.5, -sqrt(3) / 2)
    }
  );

  // An example of a function of two variables
  template<class T>
  CubicFunction<T> zCubedMinusi = CubicFunction<T>(
    // "z^3-i",
    { 
      1, 0, -3, 0, 0, 
      0, 3, 0, -1, -1
    },
    std::vector<Coords2D<T>>{
      Coords2D<T>(0, -1),
      Coords2D<T>( sqrt(3) / 2, 0.5),
      Coords2D<T>(-sqrt(3) / 2, 0.5)
    }
  );

  void fractal() {
    const int NUM_I = 1500;
    const int NUM_J = 1500;

    const int displaySize = 1;

    CubicFunction<float> f = zCubedMinus1<float>;
    const float LOW_X = -1.5;
    const float HIGH_X = 1.5;
    const float LOW_Y = -1.5;
    const float HIGH_Y = 1.5;
    const float accuracy_tolerance = 1e-2; // spectacular image!

    //CubicFunction<double> f = zCubedMinus1<double>;
    //const double LOW_X = -1.5;
    //const double HIGH_X = 1.5;
    //const double LOW_Y = -1.5;
    //const double HIGH_Y = 1.5;
    //const double accuracy_tolerance = 1e-6;

    //CubicFunction<long double> f = zCubedMinus1<long double>;
    //const long double LOW_X = -1.5;
    //const long double HIGH_X = 1.5;
    //const long double LOW_Y = -1.5;
    //const long double HIGH_Y = 1.5;
    //const long double accuracy_tolerance = 1e-6;

    //CubicFunction<long double> f = zCubedMinusi<long double>;
    //const long double LOW_X = -1.5;
    //const long double HIGH_X = 1.5;
    //const long double LOW_Y = -1.5;
    //const long double HIGH_Y = 1.5;
    //const long double accuracy_tolerance = 1e-6;

    // each known root will create a region of a different colour
    // with an additional region for non-convergent initial points
    std::vector<std::vector<Coords2D<float>>> ptsData(f.roots.size() + 1);

    assessConvergence(
      f,
      ptsData,
      NUM_I, NUM_J, 
      LOW_X, HIGH_X, LOW_Y, HIGH_Y,
      accuracy_tolerance
    );

    std::vector<std::vector<Point>> pts_to_display(f.roots.size() + 1);

    for (int i = 0; i < ptsData.size(); i++) {
      pts_to_display[i].resize(ptsData[i].size());
      auto f = [](Coords2D<float> pt) { return Point(pt.X(), pt.Y(), 0.0); };

      std::transform(ptsData[i].begin(), ptsData[i].end(), pts_to_display[i].begin(), f);

      std::cout << "pts_to_display[" << i << "].size() = " << pts_to_display[i].size() << "\n";
    }

    addPtsToView(pts_to_display, displaySize);
  }
}
