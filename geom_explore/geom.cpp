#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <tuple>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <chrono>
#include "../nr_explore/gaussj.hpp"
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
    std::optional<Vector*> first_derivative, 
    std::optional<Vector*> second_derivative
  ) const {
    Point p{
      start.X() + t * direction.X(), 
      start.Y() + t * direction.Y(), 
      start.Z() + t * direction.Z() 
    };

    if (first_derivative.has_value()) {
      first_derivative.value()->x = direction.X();
      first_derivative.value()->y = direction.Y();
      first_derivative.value()->z = direction.Z();
    }

    if (second_derivative.has_value()) {
      second_derivative.value()->x = 0;
      second_derivative.value()->y = 0;
      second_derivative.value()->z = 0;
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
    std::optional<Vector*> first_derivative, 
    std::optional<Vector*> second_derivative
  ) const {
    Point p{
      center.X() + radius * cos(t), 
      center.Y() + radius * sin(t),
      center.Z()
    };

    if (first_derivative.has_value()) {
      first_derivative.value()->x = -radius * sin(t);
      first_derivative.value()->y = radius * cos(t);
    }

    if (second_derivative.has_value()) {
      second_derivative.value()->x = -radius * cos(t);
      second_derivative.value()->y = -radius * sin(t);
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
    Point p = c.evaluate(t, &first, &second);

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
      Point p = c.evaluate(t, &firstDeriv);
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
        Point p = c.evaluate(t, &firstDeriv, &secondDeriv);
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
    Point p = c.evaluate(t, &deriv);
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
    double rootEstimate = guess;
    for (int j = 0; j < MAX_ITERATIONS; j++) {
      double f = funcd(rootEstimate);
      double df = funcd.df(rootEstimate);
      if (df == 0.0) 
        throw std::runtime_error("df == 0.0 in rtnewt");
      double dx = f / df;
      rootEstimate -= dx;
      // std::cout << "j is " << j << " and rootEstimate is " << rootEstimate << "\n";
      if ((lower_bound_param - rootEstimate) * (rootEstimate - upper_bound_param) < 0.0)
        throw std::runtime_error("Jumped out of range in rtnewt after " + std::to_string(j) + " iterations");
      if (abs(dx) < accuracy_tolerance) 
        return rootEstimate; 
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
      std::cerr << "Unknown exception caught! for f1 in test_distance_sq_fn" << std::endl;
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
        std::cerr << "Unknown exception caught! for f2 in test_distance_sq_fn" << std::endl;
      }
    }

    geom_examples::writeGeometryToJSON();
  }

  // A function of two variables uses (x,y)
  // and a function of a complex variable uses x+iy
  template <class T>
  class Coords2D {

  public:
    T x, y;

    // Constructor
    Coords2D(T x = 0, T y = 0): x(x), y(y) {}

    // Getter functions
    T X() const { return x; }
    T Y() const { return y; }

    bool operator<(const Coords2D<T>& rhs) const{

      bool lhsNan = std::isnan(X()) || std::isnan(Y());
      bool rhsNan = std::isnan(rhs.X()) || std::isnan(rhs.Y());
  
      if (lhsNan || rhsNan) {
          // If both are NaN, treat them as equal (i.e., neither is less)
          if (lhsNan && rhsNan) return false;
          // Otherwise, decide that NaN is "less" than a valid number.
          return lhsNan && !rhsNan;
      }

      return (X() < rhs.X()) || 
        ((X() == rhs.X()) && (Y() < rhs.Y()));
    }
    Coords2D<T> operator+(const Coords2D<T>& rhs) const{
      return Coords2D(x + rhs.x, y + rhs.y);
    }
    Coords2D<T> operator-(const Coords2D<T>& rhs) const{
      return Coords2D(x - rhs.x, y - rhs.y);
    }
  };

  // A certain kind of cubic function of two real values
  // inspired by the complex function z^3 - 1
  template <class T>
  class CubicFunction {
    public:
      explicit CubicFunction(
        std::vector<T> coeffs
      ) {
        //     ax^3 + bx^2y + cxy^2 + dy^3 + e
        // + i(fx^3 + gx^2y + hxy^2 + iy^3 + j)
        a = coeffs[0]; //   x^3
        b = coeffs[1]; //   x^2y
        c = coeffs[2]; //   xy^2
        d = coeffs[3]; //   y^3
        e = coeffs[4]; //   1
        f = coeffs[5]; // i x^3
        g = coeffs[6]; // i x^2y
        h = coeffs[7]; // i xy^2
        i = coeffs[8]; // i y^3
        j = coeffs[9]; // i
      }

    private: 
      // coefficients of the cubic polynomial
      T a = 0.0, b = 0.0, c = 0.0, d = 0.0;
      T e = 0.0;

      T f = 0.0, g = 0.0, h = 0.0, i = 0.0;
      T j = 0.0;

    public: 

      // expand a cubic polynomial
      void eval(
        const Coords2D<T>& z, 
        Coords2D<T>* output,
        Coords2D<T>* dfx,
        Coords2D<T>* dfy
      ) {
        // std::cout << "in cubic funtion eval, z is " << z.X() << ", " << z.Y() << "\n";
        if (output) {
          output->x = 
            a * z.X() * z.X() * z.X() // a = 1    1 * z.X() * z.X() * z.X()
          + b * z.X() * z.X() * z.Y() // b = 0    0
          + c * z.X() * z.Y() * z.Y() // c = -3  -3 * z.X() * z.Y() * z.Y()
          + d * z.Y() * z.Y() * z.Y() // d = 0    0
          + e;                        // e = -1  -1
          output->y = 
            f * z.X() * z.X() * z.X() // f = 0    0
          + g * z.X() * z.X() * z.Y() // g = 3    3 * z.X() * z.X() * z.Y()
          + h * z.X() * z.Y() * z.Y() // h = 0    0
          + i * z.Y() * z.Y() * z.Y() // i = -1   -1 * z.Y() * z.Y() * z.Y()
          + j;                        // j = 0    0
          // std::cout << "output is " << output->X() << ", " << output->Y() << "\n";
        }
        if (dfx) {
          dfx->x = 
            3 * a * z.X() * z.X() // a = 1    1 * 3 * z.X() * z.X()
          + 2 * b * z.X() * z.Y() // b = 0    0
          +     c * z.Y() * z.Y() // c = -3   -3 * z.Y() * z.Y()
          ;
          dfx->y = 
            3 * f * z.X() * z.X() // f = 0    0
          + 2 * g * z.X() * z.Y() // g = 3    2 * 3 * z.X() * z.Y()
          +     h * z.Y() * z.Y() // h = 0    0
          ;
          // std::cout << "dfx is " << dfx->X() << ", " << dfx->Y() << "\n";
        }
        if (dfy) {
          dfy->x = 
                b * z.X() * z.X() // b = 0    0
          + 2 * c * z.X() * z.Y() // c = -3   -6 * z.X() * z.Y()
          + 3 * d * z.Y() * z.Y() // d = 0    0
          ;
          dfy->y = 
                g * z.X() * z.X() // g = 3    3 * z.X() * z.X()
          + 2 * h * z.X() * z.Y() // h = 0    0
          + 3 * i * z.Y() * z.Y() // i = -1   -1 * 3 * z.Y() * z.Y()
          ;
          // std::cout << "dfy is " << dfy->X() << ", " << dfy->Y() << "\n";
        }
      }
  };


  // Every iterative root-finding process needs a way to
  // compute a "step" from one "estimate" to the next.
  // For processes with 2D input and 2D output, this is a 
  // matter of finding a 2x2 Jacobian, inverting it, and
  // multiplying it by the function value.
  template <
    class numType,    // e.g. float, double, long double
    class inputT,     // e.g. Coords2D<numType>
    class outputT,    // e.g. Coords2D<numType>
    class F           // e.g. CubicFunction<numType>
  >
  class StepFinder2D2D {

    // how to find matrix inverses
    // for 2x2 we can compute directly, for bigger matrices,
    // use the general gaussian elimination method
    const bool useGaussj;

    // internal data used to hold the data for gauss inversion
    std::vector<std::vector<numType>> a;
    std::vector<std::vector<numType>> b;

    // internal data used to hold the results of function evaluation
    outputT f, dfx, dfy;

    public: 
      StepFinder2D2D(bool useGaussj = false): useGaussj(useGaussj) {
        // set up matrices used for gaussj inversion of Jacobian
        a.resize(2);
        b.resize(2);
        a[0] = {0, 0};
        a[1] = {0, 0};
        b[0] = {0, 0};
        b[1] = {0, 0};
      }
    
      inputT findStep(
        F &funcd, 
        const inputT& rootEstimate
      ) {
        bool printDebug = false;
        funcd.eval(rootEstimate, &f, &dfx, &dfy);

        //std::cout << "f is " << f.X() << ", " << f.Y() << "\n";
        //std::cout << "dfx is " << dfx.X() << ", " << dfx.Y() << "\n";
        //std::cout << "dfy is " << dfy.X() << ", " << dfy.Y() << "\n";

        // Newton's rule
        // F is a column vector F_x, F_y
        // Jacobian J is a 2x2 matrix
        // ( d/dx F_x    d/dy F_x)
        // ( d/dx F_y    d/dy F_y)
        // 
        // x_new = x_old - J^{-1}F

        numType step_x, step_y;
        if (useGaussj) {
          b[0][0] = 1;
          b[0][1] = 0;
          b[1][0] = 0;
          b[1][1] = 1;

          a[0][0] = dfx.X();
          a[0][1] = dfx.Y();
          a[1][0] = dfy.X();
          a[1][1] = dfy.Y();

          nr_explore::gaussj(a, b);
          // nr_explore::printMatrices(a, b);
          
          step_x = (   b[0][0] * f.X() + b[1][0] * f.Y() );
          step_y = (   b[0][1] * f.X() + b[1][1] * f.Y() );
        } else {
          numType det = dfx.X() * dfy.Y() - dfx.Y() * dfy.X();
          if(det == 0) {
            throw("Zero derivative determinant in complex_newt");
          }

          step_x = (   dfy.Y() * f.X() - dfy.X() * f.Y() ) / det;
          if (printDebug) {
            std::cout << "  dfy.Y() * f.X() is   " << dfy.Y() << " * " << f.X() << "\n";
            std::cout << "- dfy.X() * f.Y() is - " << dfy.X() << " * " << f.Y() << "\n";
            std::cout << "step_x is " << step_x << "\n";
          }

          step_y = ( - dfx.Y() * f.X() + dfx.X() * f.Y() ) / det;
          if (printDebug) {
            std::cout << "- dfx.Y() * f.X() is - " << dfx.Y() << " * " << f.X() << "\n";
            std::cout << "  dfx.X() * f.Y() is   " << dfx.X() << " * " << f.Y() << "\n";
            std::cout << "step_y is " << step_y << "\n";
          }
        }

        if (printDebug) {
          std::cout << "step_x is " << step_x << "\n";
          std::cout << "step_y is " << step_y << "\n";
        }
        inputT step(step_x, step_y);
        return step;
      }
  };

  template<
    class inputT // e.g. Coords2D<double>
  >
  class RangeChecker2D {
    inputT lower_bound, upper_bound;
  public:
    RangeChecker2D(
      inputT lower_bound, 
      inputT upper_bound
    ) : lower_bound(lower_bound), upper_bound(upper_bound) {}

    bool inRange(inputT& rootEstimate) {
      bool printDebug = false;

      if (printDebug) {
        std::cout << "lower_bound is " << lower_bound.X() << ", " << lower_bound.Y() << "\n";
        std::cout << "upper_bound is " << upper_bound.X() << ", " << upper_bound.Y() << "\n";
        std::cout << "rootEstimate is " << rootEstimate.X() << ", " << rootEstimate.Y() << "\n";
        std::cout << "(lower_bound.X() - rootEstimate.X()) * (rootEstimate.X() - upper_bound.X()) = "
          << (lower_bound.X() - rootEstimate.X()) << " * " << (rootEstimate.X() - upper_bound.X()) << " = "
          << (lower_bound.X() - rootEstimate.X()) *  (rootEstimate.X() - upper_bound.X()) << "\n";
        std::cout << "(lower_bound.Y() - rootEstimate.Y()) * (rootEstimate.Y() - upper_bound.Y()) = "
          << (lower_bound.Y() - rootEstimate.Y()) << " * " << (rootEstimate.Y() - upper_bound.Y()) << " = "
          << (lower_bound.Y() - rootEstimate.Y()) *  (rootEstimate.Y() - upper_bound.Y()) << "\n";
      }
      return (lower_bound.X() - rootEstimate.X()) * (rootEstimate.X() - upper_bound.X()) > 0.0 &&
              (lower_bound.Y() - rootEstimate.Y()) * (rootEstimate.Y() - upper_bound.Y()) > 0.0;
      }
  };

  template<
    class T,
    class inputT
  >
  class ConvergenceChecker2D {
    private:
      T ptol_sq; // parameter-space tolerance-sqd
    public:
      ConvergenceChecker2D(T ptol_sq): ptol_sq(ptol_sq) {
      }
      bool hasConverged(inputT step) {
        return step.X() * step.X() + step.Y() * step.Y() < ptol_sq;
      }
  };

  // Use an iterative process to seek a root of a function
  // Templatise so that we can explore floats, doubles, long doubles
  // Inputs could be of different dimensions.
  template <
    class numType,    // e.g. float, double, long double
    class inputT,     // e.g. Coords2D<numType>
    class RangeChecker, // e.g. RangeChecker2D<numType, inputT>
    class outputT,    // e.g. Coords2D<numType>
    class F,          // e.g. CubicFunction<numType>
    class StepFinder, // e.g. StepFinder2D2D<numType, inputT, outputT, F>
    class ConvergenceChecker // e.g. ConvergenceChecker2D<numType, inputT>
  >
  inputT newtonRaphson(
    F &funcd, // needs operator(), dfx(), dfy() all returning T
    inputT& guess,
    RangeChecker& rangeChecker,
    StepFinder& stepFinder,
    ConvergenceChecker& convergenceChecker,
    const int MAX_ITERATIONS
  ) {
    const bool printDebug = false;

    // start the iteration at the provided guess
    inputT rootEstimate = guess;
    for (int j = 0; j < MAX_ITERATIONS; j++) {

      // use the provided StepFinder to work out where to go next
      inputT step = stepFinder.findStep(funcd, rootEstimate);
      if (printDebug) {
        std::cout << "step is " << step.X() << ", " << step.Y() << "\n";
      }
      // apply the step
      rootEstimate = rootEstimate - step;
      if (printDebug) {
        std::cout << "j is " << j << " and rootEstimate is " << rootEstimate.X() << ", " << rootEstimate.Y() << "\n";
      }

      // use the provided RangeChecker to see if we have jumped out of range
      // and should stop
      if (!rangeChecker.inRange(rootEstimate)) {
        if (printDebug) {
          std::cout << "Out of range\n";
        }
        throw("Jumped out of range in newtonRaphson");
      }

      // use the provided ConvergenceChecker to see if we have converged
      // and should stop
      if (convergenceChecker.hasConverged(step)) {
        if (printDebug) {
          std::cout << "Converged to " << rootEstimate.X() << ", " << rootEstimate.Y() << " after " << j << " iterations\n";
        }
        return rootEstimate; 
      }
    }
    throw("Maximum number of iterations exceeded in newtonRaphson");
  }

  template <class T>
  class ColorPatch2D {
      
  protected:
    Coords2D<T> solution;

  public:
    uint32_t color;

    // Constructor
    ColorPatch2D(
      Coords2D<T> sol,
      uint32_t col
    ): solution(sol), color(col) {}

    T X() const { return solution.X(); }
    T Y() const { return solution.Y(); }

    bool operator<(const ColorPatch2D<T>& rhs) const{
      return solution < rhs.solution;
    }
  };

  auto colors = generateMutedColors(2000);

  // Add the start to the foundSolutions map, e.g.
  // convergence to a value within tolerance of a value already in the map
  //   adds the start to the vector which is the map value
  // convergence to a value not within tolerance of a value already in the map
  //   adds a new key value pair with the start in the value
  // no convergence 
  //   adds the start to the vector which is the map value for a key of NaNs
  
  template <class T>
  void addToMap(
    Coords2D<T>& guess,
    Coords2D<T>& newtonResult,
    std::function<bool(ColorPatch2D<T>, Coords2D<T>)>& patchMatcher,
    std::map<ColorPatch2D<T>, std::vector<Coords2D<T>>>& foundSolutions
  ){
    bool printDebug = false;

    bool addedToMap = false;
    for (auto& kv : foundSolutions) {
      auto& key = kv.first;
      if (printDebug) {
        std::cout << "compare with (" << key.X() << ", " << key.Y() << ")\n";
      }

      bool matchedRoot = patchMatcher(key, newtonResult);

      if (matchedRoot) {
        if (printDebug) {
          std::cout << "add to existing collection for key (" 
            << key.X() << ", " << key.Y() << ")\n";
        }
  
        kv.second.push_back(guess);
        addedToMap = true;
        break;
      }
    }
  }

  // Use 2D Newton Raphson for a given function f
  // starting from a given guess.
  template <class T>
  Coords2D<T> NR2DConverges(
    CubicFunction<T>& f,
    Coords2D<T>& guess,
    T tol_sq,  // parameter-space tol-sqd
    StepFinder2D2D<T, Coords2D<T>, Coords2D<T>, CubicFunction<T>>& stepFinder
  ) {
    bool printDebug = false;

    Coords2D<T> newtonResult(std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN());
    try {
      using InputType = Coords2D<T>;
      using RangeCheckerType = RangeChecker2D<InputType>;
      using OutputType = Coords2D<T>;
      using FunctionType = CubicFunction<T>;
      using ConvergenceType = ConvergenceChecker2D<T, InputType>;

      InputType start = guess;

      RangeCheckerType rangeChecker(
        InputType(-100.0, -100.0),
        InputType(100.0, 100.0) 
      );
 
      ConvergenceType convergence(tol_sq);  // parameter-space tol-sqd

      //std::cout << "newtonResult on f from guess " << start.X() << " + i * " << start.Y() << "\n";
      newtonResult = newtonRaphson<
        T,
        InputType,
        RangeCheckerType,
        OutputType,
        FunctionType,
        StepFinder2D2D<T, Coords2D<T>, Coords2D<T>, CubicFunction<T>>,
        ConvergenceType
      >(
        f,
        start,
        rangeChecker,
        stepFinder,
        convergence,
        100
      );
      // std::cout << newtonResult.X() << " + i * " << newtonResult.Y() << "\n";

    } catch (const char* msg) {
      // std::cerr << "Error: " << msg << std::endl;
    } catch (const std::exception& e) {
      std::cerr << "Exception caught in NR2DConverges: " << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Unknown exception caught! in NR2DConverges" << std::endl;
    }

    if (printDebug) {
      std::cout << "start (" << guess.X() << ", " << guess.Y()
        << ") converged to (" << 
        newtonResult.X() << ", " << newtonResult.Y() << ")\n";
    }

    // std::cout << "foundSolutions.size() = " << foundSolutions.size() << "\n";
    return newtonResult;
  }

  template<class T>
  bool coordsMatch(
    const Coords2D<T>& lhs, 
    const Coords2D<T>& rhs,
    T tol
  ) {
    return (
      std::isnan(lhs.X()) && std::isnan(rhs.X())
    ) || (abs(lhs.X() - rhs.X()) < tol
       && abs(lhs.Y() - rhs.Y()) < tol);  
  }

  template<class T>
  bool patchMatchesCoords(
    const ColorPatch2D<T>& lhs, 
    const Coords2D<T>& rhs,
    T tol_sq
  ) {
    if (std::isnan(lhs.X()) && std::isnan(rhs.X())) {
      return true;
    }
    const T x_step = lhs.X() - rhs.X();
    const T y_step = lhs.Y() - rhs.Y();
    return (x_step * x_step + y_step * y_step) < tol_sq;
  }

  template<class T>
  bool patchMatchesCoordsFunky(
    const ColorPatch2D<T>& lhs, 
    const Coords2D<T>& rhs,
    T tol
  ) {
    const bool printDebug = false;

    // This is a pretty bizarre way to ask
    // "Have we converged to this key"
    // ie. this known result?
    // By setting a very strict tolerance
    // for the y-coordinate
    // some points say "no match" and get omitted from the plot
    const bool matched = (
      std::isnan(lhs.X()) && std::isnan(rhs.X())
    ) || (abs(lhs.X() - rhs.X()) < tol
      && abs(lhs.Y() - rhs.Y()) < 1e-6); 
    if (printDebug) {
      std::cout << "lhs (" << lhs.X() << ", " << lhs.Y() << ") "
        << "rhs (" << rhs.X() << ", " << rhs.Y() << ") "
        << "matched is " << matched << "\n";
    }
    return matched;  
  }

  template<class T>
  bool patchMatchesCoordsFunky2(
    const ColorPatch2D<T>& lhs, 
    const Coords2D<T>& rhs,
    T tol
  ) {
    // This is a pretty bizarre way to ask
    // "Have we converged to this key"
    // ie. this known result?
    // By setting a very strict tolerance
    // for the y-coordinate
    // some points say "no match" and get omitted from the plot
    if (std::isnan(lhs.X()) && std::isnan(rhs.X())) {
      return true;
    }
    const T xStep = lhs.X() - rhs.X();
    const T yStep = lhs.Y() - rhs.Y();
    return xStep * xStep + yStep * yStep < 1e-11;
  }

  template<class T>
  bool patchesMatch(
    const ColorPatch2D<T>& lhs, 
    const ColorPatch2D<T>& rhs, 
    T tol
  ) {
    return (
      std::isnan(lhs.X()) && std::isnan(rhs.X())
    ) || (abs(lhs.X() - rhs.X()) < tol
       && abs(lhs.Y() - rhs.Y()) < tol);  
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
    std::map<ColorPatch2D<T>, std::vector<Coords2D<T>>>& foundSolutions,
    std::function<bool(ColorPatch2D<T>, Coords2D<T>)>& patchMatcher,
    int NUM_I, 
    int NUM_J, 
    T LOW_X, 
    T HIGH_X, 
    T LOW_Y, 
    T HIGH_Y,
    T tol_sq  // parameter-space tol-sqd
  ) {
    bool printDebug = false;
    std::mutex mtx; // protects access to ptsConverged and ptsNotConverged

    // Lambda that processes a block of rows (i values)
    auto process_range = [&](int start_i, int end_i) {
      std::map<ColorPatch2D<T>, std::vector<Coords2D<T>>> localFoundSolutions = foundSolutions;

      // Initialise one StepFinder and internal Coords data
      // for use in all the following iterations
      StepFinder2D2D<T, Coords2D<T>, Coords2D<T>, CubicFunction<T>> stepFinder;
      for (int i = start_i; i < end_i; i++) {
        for (int j = 0; j < NUM_J; j++) {
          Coords2D<T> start(
            LOW_X + (HIGH_X - LOW_X) / NUM_I * i,
            LOW_Y + (HIGH_Y - LOW_Y) / NUM_J * j
          );
          Coords2D<T> newtonResult = NR2DConverges<T>(
            f, 
            start,
            tol_sq,  // parameter-space tol-sqd
            stepFinder
          );
          addToMap(
            start,
            newtonResult, 
            patchMatcher,
            localFoundSolutions
          );
        }
      }
      // Merge local results into the global map with locking
      std::lock_guard<std::mutex> lock(mtx);
      for (const auto& localKv : localFoundSolutions) {
        auto& localKey = localKv.first;

        // look for a matching key value in foundSolutions
        bool addedToMap = false;
        for (auto& foundKv : foundSolutions) {
          auto& foundKey = foundKv.first;
          if (patchesMatch(localKey, foundKey, tol_sq)) {
            if (printDebug) {
              std::cout << "add to existing collection for key (" << foundKey.X() << ", " << foundKey.Y() << ")\n";
            }

            // add the contents of localKv.second to foundKey.second
            foundKv.second.insert(foundKv.second.end(), localKv.second.begin(), localKv.second.end());
            addedToMap = true;
            break;
          }
        }
        if (!addedToMap) {
          // add a new key value pair to the map
          //std::cout << "create new collection for newtonResult (" << localKey.X() << ", " << localKey.Y() << ")\n";
          foundSolutions[localKey] = localKv.second;
        }        
      }
    };

    // Determine the number of threads available.
    unsigned int numThreads = std::thread::hardware_concurrency();
    if (numThreads == 0) { numThreads = 2; } // Fallback in case hardware_concurrency returns 0

    if (printDebug) {
      numThreads = 1;
      std::cout << "Using " << numThreads << " threads\n";
    }

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

  template<class T>
  void addPtsToView(
    std::map<ColorPatch2D<T>, std::vector<Coords2D<T>>>& foundSolutions,
    int displaySize
  ){
    std::vector<PtCollection> ptsColls;

    for (auto& kv : foundSolutions) {
      std::vector<Point> pts;
      const std::vector<Coords2D<T>>& starts = kv.second;
      for (const auto& pt : starts) {
        pts.push_back(Point(pt.X(), pt.Y(), 0.0));
      }
      // std::cout << "pts_to_display elt .size() = " << pts.size() << "\n";

      PtCollection ptsObj = {
        .displaySize = displaySize,
        .color = kv.first.color,
        .pts = pts,
        .isLine = false
      };
      ptsColls.push_back(ptsObj);
    }

    // std::cout << "adding " << ptsColls.size() << " collection(s) to view\n";
    addGeometryToView(ptsColls);    
  }

  // An example of a function of two variables
  template<class T>
  CubicFunction<T> zCubedMinus1 = CubicFunction<T>(
    // "z^3-1",
    { 
      1, 0, -3, 0, -1, 
      0, 3, 0, -1, 0
    }
  );
  
  class CurveDifference2D {
    private:
      Curve& c1;
      Curve& c2;
    public:
      using InputType = Coords2D<double>; // think parameter-pair (s, t)
      using OutputType = Coords2D<double>; // think space position (x, y)
      explicit CurveDifference2D(Curve& c1, Curve& c2)
        : c1(c1), c2(c2) 
        { }

      void eval(
        const Coords2D<double>& z, 
        Coords2D<double>* output,
        Coords2D<double>* dfx,
        Coords2D<double>* dfy
      ) {
        if (output) {
          Point pos1 = c1.evaluate(z.X());
          Point pos2 = c2.evaluate(z.Y()); 
          output->x = pos1.X() - pos2.X();
          output->y = pos1.Y() - pos2.Y();
        }
        if (dfx) {
          Vector first1;
          c1.evaluate(z.X(), &first1);
          dfx->x = first1.X();
          dfx->y = first1.Y();
        }
        if (dfy) {
          Vector first2;
          c2.evaluate(z.Y(), &first2);
          dfy->x = -first2.X();
          dfy->y = -first2.Y();
        }
      }
  };

  void testNrZCubedMinus1() {
    // function of a complex variable f(z) = z^3 - 1
    // expect convergence if closeish to one of three roots, other convergence
    // from other starting points behaves in a fractal-like way

    using InputType = Coords2D<float>;
    using RangeCheckerType = RangeChecker2D<InputType>;
    using OutputType = Coords2D<float>;
    using FunctionType = CubicFunction<float>;
    using StepFinderType = StepFinder2D2D<float, InputType, OutputType, FunctionType>;
    using ConvergenceType = ConvergenceChecker2D<float, InputType>;

    InputType startZCubedMinus1 = InputType(0, -1);

    RangeCheckerType rangeChecker(
      InputType(-100.0, -100.0),
      InputType(100.0, 100.0) 
    );

    StepFinderType stepFinder;

    ConvergenceType convergence(1e-3);

    Coords2D<float> newtonResultzCubedMinus1 = newtonRaphson<
      float,
      InputType,
      RangeCheckerType,
      OutputType,
      FunctionType,
      StepFinderType,
      ConvergenceType
  >(
      zCubedMinus1<float>, 
      startZCubedMinus1,
      rangeChecker,
      stepFinder,
      convergence,  // parameter-space tol-sqd
      100
    );

    std::cout << "z^3 - 1 starting from " <<  startZCubedMinus1.X() << ", " << startZCubedMinus1.Y() << " yielded " 
      << newtonResultzCubedMinus1.X() << ", " << newtonResultzCubedMinus1.Y() << "\n";
  }

  Point intersect2DCurves(
    Curve& c1, 
    Curve& c2,
    Coords2D<double>& start
  ) {

    using InputType = Coords2D<double>;
    using RangeCheckerType = RangeChecker2D<InputType>;
    using OutputType = Coords2D<double>;
    using FunctionType = CurveDifference2D;
    using ConvergenceType = ConvergenceChecker2D<double, InputType>;


    RangeCheckerType rangeChecker(
      InputType(-100.0, -100.0),
      InputType(100.0, 100.0)
    );

    StepFinder2D2D<double, InputType, OutputType, FunctionType> stepFinder;

    ConvergenceType convergence(1e-3);

    // function representing the difference between two lines
    // expect (immediate) convergence to their intersection
    CurveDifference2D diff(c1, c2);



    //Double2DCoords startAxes = Double2DCoords(0, -1);
    InputType newtonResultl1l2 = newtonRaphson<
      double,
      InputType,
      RangeCheckerType,
      OutputType,
      FunctionType,
      StepFinder2D2D<double, InputType, OutputType, FunctionType>,
      ConvergenceType
    >(
      diff, 
      start,
      rangeChecker,
      stepFinder,
      convergence,  // parameter-space tol-sqd
      100
    );
    std::cout << "iteration starting from " <<  start.X() << ", " << start.Y() << " yielded " 
      << newtonResultl1l2.X() << ", " << newtonResultl1l2.Y() << "\n";

    Point onC1 = c1.evaluate(newtonResultl1l2.X());
    Point onC2 = c2.evaluate(newtonResultl1l2.Y());

    std::cout << "onL1 = " <<  onC1.X() << ", " << onC1.Y() << "\n";
    std::cout << "onL2 = " <<  onC2.X() << ", " << onC2.Y() << "\n";

    return onC1; // TOO MUCH TO DO - what if we want to return onC2? Convergence fails?
  }

  void intersectAxes2D() {
    Line l1(Point(1,2,0), Vector(1,0,0));
    Line l2(Point(22,11,0), Vector(0,1,0));
    Coords2D<double> start = Coords2D<double>(1, 1);
    intersect2DCurves(
      l1,
      l2,
      start
    );
  }

  void intersectLines2D() {
    Line l1(Point(1,2,0), Vector(1,1,0));
    Line l2(Point(22,11,0), Vector(2,1,0));
    Coords2D<double> start = Coords2D<double>(0, -1);
    intersect2DCurves(
      l1,
      l2,
      start
    );
  }

  void intersectCircles2D() {

    Circle c1(Point(0.0, 0.0, 0.0), 2);
    Circle c2(Point(1.0, 2.0, 0.0), sqrt(5.0));
    Coords2D<double> start = Coords2D<double>(0, -1);

    addCurveToView(c1, 0.0, 2*M_PI, GREEN, 2);
    addCurveToView(c2, 0.0, 2*M_PI, RED, 2);

    Point onC1 = intersect2DCurves(
      c1,
      c2,
      start
    );

    addPointToView(onC1, YELLOW, 7);
  }

  void testNewtonRaphson2Dinput(){

    try {
      testNrZCubedMinus1();

      Line xaxis(Point(0,0,0), Vector(1,0,0));
      Line yaxis(Point(0,0,0), Vector(0,1,0));
      Line zaxis(Point(0,0,0), Vector(0,0,1));
  
      addCurveToView(xaxis, 0.0, 1.0, RED, 2);
      addCurveToView(yaxis, 0.0, 1.0, GREEN, 2);
      addCurveToView(zaxis, 0.0, 1.0, BLUE, 2);

      intersectAxes2D();    // should find (0,0)
      intersectLines2D();   // should converge immediately
      intersectCircles2D(); // a more intereseting case

    } catch (const std::runtime_error& e) {
      std::cerr << "in testNewtonRaphson2Dinput got an error: " << e.what() << std::endl;
    } catch (const char* msg) {
      std::cerr << "in testNewtonRaphson2Dinput got an error: " << msg << std::endl;
    } catch (...) {
      std::cerr << "in testNewtonRaphson2Dinput got an Unknown exception caught!" << std::endl;
    }

  }

  template <class T>
  void findSolutions(
    CubicFunction<T>& f,
    T LOW_X,
    T HIGH_X,
    T LOW_Y,
    T HIGH_Y,
    T accuracy_tolerance, 
    std::map<ColorPatch2D<T>, std::vector<Coords2D<T>>>& foundSolutions
  ) {
    bool printDebug = false;
    std::vector<Coords2D<T>> solutions = {};

    StepFinder2D2D<T, Coords2D<T>, Coords2D<T>, CubicFunction<T>> stepFinder;

    // find a number of accurate solutions for f
    const int NUM_SAMPLES = 10;
    for (int i = 0; i < NUM_SAMPLES; i++) {
      for (int j = 0; j < NUM_SAMPLES; j++) {
        Coords2D<T> start(
          LOW_X + (HIGH_X - LOW_X) / NUM_SAMPLES * i,
          LOW_Y + (HIGH_Y - LOW_Y) / NUM_SAMPLES * j
        );
        Coords2D<T> newtonResult = NR2DConverges<T>(
          f, 
          start,
          1e-6,
          stepFinder
        );
        /*
        std::cout << "solution = (" 
          << newtonResult.X() << ", " << newtonResult.Y() 
          << ") from (" 
          << start.X() << ", " << start.Y() 
          << ")\n";
        */
        bool matchedSol = false;
        for(auto& knownSol: solutions) {
          if (coordsMatch(newtonResult, knownSol, accuracy_tolerance)) {
            matchedSol = true;
            break;
          }
        }
        if (!matchedSol) {
          solutions.push_back(newtonResult);
        }
      }
    }

    for (int i = 0; i < solutions.size(); i++) {
      const auto& knownSol = solutions[i];
      if (std::isnan(knownSol.X())) {
        foundSolutions[ColorPatch2D<T>(
          knownSol,
          WHITE
        )] = {};  
        continue;
      }
      if (i >= colors.size()) {
        //if (printDebug) {
          std::cout << "not plotted solution = (" << knownSol.X() << ", " << knownSol.Y() << ")\n";
        //}
        continue;    
      }
      const auto& col = colors[i];
      //if (printDebug) {
        std::cout << "solution = (" << knownSol.X() << ", " << knownSol.Y() << ") plotted "<< col << " \n";
      //}
      foundSolutions[ColorPatch2D<T>(
        knownSol,
        col
      )] = {};
    }
  }

  template <class T> 
  void makePlot(
    CubicFunction<T>& f,
    std::map<ColorPatch2D<T>, std::vector<Coords2D<T>>>& foundSolutions,
    T LOW_X,
    T HIGH_X,
    T LOW_Y,
    T HIGH_Y,
    T accuracy_tolerance,
    std::function<bool(ColorPatch2D<T>, Coords2D<T>)>& patchMatcher,
    int NUM_I,
    int NUM_J,
    int displaySize
  ) {

    // Start timer
    auto start = std::chrono::high_resolution_clock::now();

    if (foundSolutions.size() == 0) {
      findSolutions(
        f,
        LOW_X,
        HIGH_X,
        LOW_Y,
        HIGH_Y,
        accuracy_tolerance, 
        foundSolutions
      );  
    }

    assessConvergence(
      f,
      foundSolutions,
      patchMatcher,
      NUM_I, NUM_J, 
      LOW_X, HIGH_X, LOW_Y, HIGH_Y,
      accuracy_tolerance * accuracy_tolerance
    );

    // End timer
    auto end = std::chrono::high_resolution_clock::now();
    // Calculate duration
    std::chrono::duration<double> duration = end - start;
    std::cout << "time for fractal " << duration.count() << "\n";

    struct PBMPt {
      double x;
      double y;
      uint32_t color;
    };
  
    std::vector<PBMPt> pbmPts;
    for (auto& kv : foundSolutions) {
      const auto& key = kv.first;
      const auto& pts = kv.second;
      for (const auto& pt : pts) {
        pbmPts.push_back(
          { 
            .x = pt.X(), 
            .y = pt.Y(), 
            .color = key.color
          }
        );
      }
    }
    // sort pbmPts according to y coord then x coord
    std::sort(pbmPts.begin(), pbmPts.end(), [](const auto& lhs, const auto& rhs) {
      if (lhs.y == rhs.y) {
        return lhs.x < rhs.x;
      }
      return lhs.y < rhs.y;
    });
    // open output file for writing
    std::ofstream outputFile("viewer/public/output/fractal.pbm");
    if (!outputFile) {
      std::cerr << "Error opening output file for writing\n";
      return;
    }
    // write PBM header
    outputFile << "P3\n";
    // use NUM_I and NUM_J to set the size of the image
    outputFile << NUM_I << " " << NUM_J << "\n";
    outputFile << "255\n";
    // write pixel data
    for (const auto& pt : pbmPts) {
      // turn the uint32_t representation of color
      // into a 3-tuple of RGB values
      uint32_t color = pt.color;
      uint32_t r = (color >> 16) & 0xFF;
      uint32_t g = (color >> 8) & 0xFF;
      uint32_t b = color & 0xFF;
      // write the pixel color
      outputFile << r << " " << g << " " << b << "\n";
    }

    clearView();
    addPtsToView(foundSolutions, displaySize);
  }

  template <class T>
  void fractalTyped(double origin_r, double origin_i, double radius, int num_r, int num_i) {
    std::cout << "start timing for fractal\n";

    CubicFunction<T> f = zCubedMinus1<T>;
    std::map<ColorPatch2D<T>, std::vector<Coords2D<T>>> foundSolutions;
    int displaySize = 5;

    try {
      int NUM_I = 300;
      int NUM_J = 30;
      T accuracy_tolerance = 0.001;
      T LOW_X = -10.0;
      T HIGH_X = 10.0;
      T LOW_Y = -10.0;
      T HIGH_Y = 10.0;
       // parameter-space tol-sqd
      T tol_sq = accuracy_tolerance * accuracy_tolerance;

      std::function<bool(ColorPatch2D<T>, Coords2D<T>)> patchMatcher =
        [&tol_sq](ColorPatch2D<T> patch, Coords2D<T> pos) -> bool {
          return patchMatchesCoords(patch, pos, tol_sq);
      };

      const int exampleNumber = 0;

      if (exampleNumber == 0) {

        // how many points to sample in each direction
        NUM_I = num_r; // real
        NUM_J = num_i; // imaginary
        // size of the rectangle plotted for each sample point
        displaySize = 1;

        // parameter range of the (x, y) initial starting point for iteration
        LOW_X  = origin_r - radius;
        HIGH_X = origin_r + radius;
        LOW_Y  = origin_i - radius;
        HIGH_Y = origin_i + radius;

        // what constitutes convergence to a known solution
        accuracy_tolerance = 0.012; // spectacular image!

        // a "patch" is a set of initial values that converge to the same solution
        patchMatcher =
          [&accuracy_tolerance](ColorPatch2D<T> patch, Coords2D<T> pos) -> bool {
            return patchMatchesCoordsFunky(patch, pos, accuracy_tolerance);
        };

        // what function are we iterating as z = x + i * y -> f(z)
        f = zCubedMinus1<T>;

        // prime foundSolutions with known roots if you know them
        foundSolutions[ColorPatch2D<T>(
          Coords2D<T>(1, 0),
          RED
        )] = {};
        foundSolutions[ColorPatch2D<T>(
          Coords2D<T>(-0.5, sqrt(3) / 2),
          YELLOW
        )] = {};
        foundSolutions[ColorPatch2D<T>(
          Coords2D<T>(-0.5, -sqrt(3) / 2),
          GREEN
        )] = {};

        // this does the work of scanning across the parameter space
        // and finding the patches for our solutions
        makePlot(
          f,
          foundSolutions,
          LOW_X,
          HIGH_X,
          LOW_Y,
          HIGH_Y,
          accuracy_tolerance,
          patchMatcher,
          NUM_I,
          NUM_J,
          displaySize
        );
      } else if(exampleNumber == 1) {

        NUM_I = 1000;
        NUM_J = 1000;
        displaySize = 1;

        LOW_X = -2.0;
        HIGH_X = 4.5;
        LOW_Y = -4.5;
        HIGH_Y = 2.0;

        accuracy_tolerance = 0.024;

        patchMatcher =
          [&accuracy_tolerance](ColorPatch2D<T> patch, Coords2D<T> pos) -> bool {
            return patchMatchesCoordsFunky(patch, pos, accuracy_tolerance);
        };

        f = CubicFunction<T>(
          { 
            1, 1, 1, 1, 4,
            0, 3, 0, -1, 1
          }
        );

        makePlot(
          f,
          foundSolutions,
          LOW_X,
          HIGH_X,
          LOW_Y,
          HIGH_Y,
          accuracy_tolerance,
          patchMatcher,
          NUM_I,
          NUM_J,
          displaySize
        );

      } else if(exampleNumber == 2) {

        NUM_I = 1000;
        NUM_J = 1000;
        displaySize = 1;

        LOW_X = -2.0;
        HIGH_X = 4.5;
        LOW_Y = -4.5;
        HIGH_Y = 2.0;

        accuracy_tolerance = 0.024;

        patchMatcher =
          [&accuracy_tolerance](ColorPatch2D<T> patch, Coords2D<T> pos) -> bool {
            return patchMatchesCoordsFunky2(patch, pos, accuracy_tolerance);
        };

        f = CubicFunction<T>(
          { 
            1, 1, 1, 1, 4,
            0, 3, 0, -1, 1
          }
        );

        makePlot(
          f,
          foundSolutions,
          LOW_X,
          HIGH_X,
          LOW_Y,
          HIGH_Y,
          accuracy_tolerance,
          patchMatcher,
          NUM_I,
          NUM_J,
          displaySize
        );
      } else if(exampleNumber == 3) {

        NUM_I = 1000;
        NUM_J = 1000;
        displaySize = 1;

        LOW_X = -2.0;
        HIGH_X = 4.5;
        LOW_Y = -4.5;
        HIGH_Y = 2.0;

        // when the step gets smaller than this,
        // Newton Raphson says we have converged
        accuracy_tolerance = 0.015;

        patchMatcher = [](ColorPatch2D<T> patch, Coords2D<T> pos) -> bool {
          if (std::isnan(patch.X()) && std::isnan(pos.X())) {
            return true;
          }
          const T xStep = patch.X() - pos.X();
          const T yStep = patch.Y() - pos.Y();

          // sols are expected to be within accuracy_tolerance
          // but if they're not within a tighter tol, 
          // reject them i.e. color them black

          return xStep * xStep + yStep * yStep < 1e-11; 
        };

        f = CubicFunction<T>(
          { 
            1, 1, 1, 1, 4,
            0, 3, 0, -1, 1
          }
        );
        makePlot(
          f,
          foundSolutions,
          LOW_X,
          HIGH_X,
          LOW_Y,
          HIGH_Y,
          accuracy_tolerance,
          patchMatcher,
          NUM_I,
          NUM_J,
          displaySize
        );
      } else if(exampleNumber == 4) {

        NUM_I = 1000;
        NUM_J = 1000;
        displaySize = 1;

        LOW_X = -0.5;
        HIGH_X = 1.5;
        LOW_Y = -1;
        HIGH_Y = 1;

        // when the step gets smaller than this,
        // Newton Raphson says we have converged
        accuracy_tolerance = 0.012;

        patchMatcher = [](ColorPatch2D<T> patch, Coords2D<T> pos) -> bool {
          if (std::isnan(patch.X()) && std::isnan(pos.X())) {
            return true;
          }
          const T xStep = patch.X() - pos.X();
          const T yStep = patch.Y() - pos.Y();

          // sols are expected to be within accuracy_tolerance
          // but if they're not within a tighter tol, 
          // reject them i.e. color them black

          return xStep * xStep + yStep * yStep < 1e-12; 
        };

        // this loop is aiming at adjusting the task to produce a
        // sequence of images we can animate
        for (T a = 0.0; a < 0.5; a += 0.1) {
          std::cout << "a = " << a << "\n";

          foundSolutions.clear();
          f = CubicFunction<T>(
            { 
              a,  0,   a, 0, 1,
              0,  -1,  0, 1, 0
            }
          );

          makePlot(
            f,
            foundSolutions,
            LOW_X,
            HIGH_X,
            LOW_Y,
            HIGH_Y,
            accuracy_tolerance,
            patchMatcher,
            NUM_I,
            NUM_J,
            displaySize
          );
        }
      }

    } catch (...) {

    }

  }

    void fractal(double origin_r, double origin_i, double radius, int num_r, int num_i) {
      fractalTyped<double>(origin_r, origin_i, radius, num_r, num_i);
  }

  // A function of three variables uses (x,y, z)
  // and a function of a complex variable uses x+iy
  template <class T>
  class Coords3D {
    T x, y, z;

    public:
    // Constructor
    Coords3D(T x = 0, T y = 0, T z = 0);

    // Getter functions
    T X() const;
    T Y() const;
    T Z() const;
  };

  template <class T>
  Coords3D<T>::Coords3D(T x, T y, T z) : x(x), y(y), z(z) {}

  template <class T>
  T Coords3D<T>::X() const { return x; }
  template <class T>
  T Coords3D<T>::Y() const { return y; }
  template <class T>
  T Coords3D<T>::Z() const { return z; }

  template <class T>
  Coords3D<T> operator+(const Coords3D<T>& a, const Coords3D<T>& b) {
    return Coords3D(a.x + b.x, a.y + b.y, a.z + b.z);
  }

  template <class T>
  Coords3D<T> operator-(const Coords3D<T>& a, const Coords3D<T>& b) {
    return Coords3D(a.x - b.x, a.y - b.y, a.z - b.z);
  }

  // Define operator< so that Coords3D objects can be ordered.
  template <class T>
  bool operator<(const Coords3D<T>& lhs, const Coords3D<T>& rhs) {
    bool lhsNan = std::isnan(lhs.X()) || std::isnan(lhs.Y()) || std::isnan(lhs.Z());
    bool rhsNan = std::isnan(rhs.X()) || std::isnan(rhs.Y()) || std::isnan(rhs.Z());

    if (lhsNan || rhsNan) {
        // If both are NaN, treat them as equal (i.e., neither is less)
        if (lhsNan && rhsNan) return false;
        // Otherwise, decide that NaN is "less" than a valid number.
        return lhsNan && !rhsNan;
    }

    return (lhs.X() < rhs.X()) 
        || ((lhs.X() == rhs.X()) 
          && (lhs.Y() < rhs.Y()))
        || ((lhs.X() == rhs.X()) 
          && (lhs.Y() == rhs.Y()) 
          && (lhs.Z() < rhs.Z()));
  }

  class CurveDifference3D {
    private:
      Curve& c1;
      Curve& c2;
    public:
      using InputType = Coords2D<double>; // think parameter-pair (s, t)
      using OutputType = Coords3D<double>; // think space position (x, y, z)
      explicit CurveDifference3D(Curve& c1, Curve& c2)
        : c1(c1), c2(c2) 
        { }

        OutputType operator() (const InputType& z) {
        Point pos1 = c1.evaluate(z.X());
        Point pos2 = c2.evaluate(z.Y()); 

        //std::cout << "pos1 is " << pos1.X() << " " << pos1.Y() << " " << pos1.Z() << "\n";
        //std::cout << "pos2 is " << pos2.X() << " " << pos2.Y() << " " << pos2.Z() << "\n";

        //addPointToView(pos1, RED, 7);
        //addPointToView(pos2, GREEN, 7);

        return OutputType(
          pos1.X() - pos2.X(),
          pos1.Y() - pos2.Y(),
          pos1.Z() - pos2.Z()
        );
      }
      OutputType dfx(const InputType& z) {
        Vector first1;
        c1.evaluate(z.X(), &first1);

        return OutputType(
          first1.X(),
          first1.Y(),
          first1.Z()
        );
      }
      OutputType dfy(const InputType& z) {
        Vector first2;
        c2.evaluate(z.Y(), &first2);

        return OutputType(
          -first2.X(),
          -first2.Y(),
          -first2.Z()
        );
      }
  };

  // Every iterative root-finding process needs a way to
  // compute a "step" from one "estimate" to the next.
  // For processes with 2D input and 3D output, we find a rectangular
  // Jacobian, a transpose of it, and an inverse.

  template <
    class numType,    // e.g. float, double, long double
    class inputT,     // e.g. Coords3D<numType>
    class outputT,    // e.g. Coords3D<numType>
    class F           // e.g. CurveDiff3d<numType>
  >
  class StepFinder2D3D {
    std::vector<std::vector<numType>> a;
    std::vector<std::vector<numType>> b;

    public: 
      StepFinder2D3D() {
        // set up matrices used for gaussj inversion of Jacobian
        a.resize(2);
        b.resize(2);
        a[0] = {0, 0};
        a[1] = {0, 0};
        b[0] = {0, 0};
        b[1] = {0, 0};
      }

      inputT findStep(
        F &funcd, 
        const inputT& rootEstimate
      ) {
        bool printDebug = false;
        outputT f = funcd(rootEstimate);
        outputT dfx = funcd.dfx(rootEstimate);
        outputT dfy = funcd.dfy(rootEstimate);

        if (printDebug) {
          std::cout << std::fixed << std::setprecision(6) << std::setw(5);

          std::cout << "f is \n";
          std::cout << f.X() << "\t" << f.Y() << "\t" << f.Z() << "\n";

          std::cout << "dfx is \n";
          std::cout << dfx.X() << "\t" << dfx.Y() << "\t" << dfx.Z() << "\n";

          std::cout << "dfy is \n";
          std::cout << dfy.X() << "\t" << dfy.Y() << "\t" << dfy.Z() << "\n";

          std::cout << std::defaultfloat;
        }


        // Newton's rule
        //
        // F is a 3x1 matrix; a column vector F_x, F_y, F_z
        //
        // Jacobian J is a 3x2 matrix
        // ( d/dx F_x    d/dy F_x )
        // ( d/dx F_y    d/dy F_y )
        // ( d/dx F_z    d/dy F_z )
        // 
        // We'll apply
        // x_new = x_old - (JtJ)^{-1}JtF
        //
        // Jt is a 2x3 matrix
        // ( d/dx F_x    d/dx F_y   d/dx F_z )
        // ( d/dy F_x    d/dy F_y   d/dy F_z )
        //
        // JtJ is a 2x2 matrix - call it a for inversion

        a[0][0] = dfx.X() * dfx.X() + dfx.Y() * dfx.Y() + dfx.Z() * dfx.Z();
        a[0][1] = dfx.X() * dfy.X() + dfx.Y() * dfy.Y() + dfx.Z() * dfy.Z();
        a[1][0] = dfy.X() * dfx.X() + dfy.Y() * dfx.Y() + dfy.Z() * dfx.Z();
        a[1][1] = dfy.X() * dfy.X() + dfy.Y() * dfy.Y() + dfy.Z() * dfy.Z();

        // 
        // (JtJ)^{-1} is a 2x2 matrix
        //
        // (JtJ)^{-1}Jt is 2x3 matrix
        //
        // (JtJ)^{-1}Jt F is 2x1 matrix

        // J^{-1} = 1/(det J) * ( d/dy F_y    -d/dy F_x)
        //                      (-d/dx F_y     d/dx F_x)
        // J^{-1}F = 1/(det J) * ( d/dy F_y * F_x - d/dy F_x * F_y)
        //                       (-d/dx F_y * F_x + d/dx F_x * F_y)
        // J^{-1}F = 1/(det J) * ( d/dy F_y * F_x - d/dy F_x * F_y)
        //                       (-d/dx F_y * F_x + d/dx F_x * F_y)

        b[0][0] = 1;
        b[0][1] = 0;
        b[1][0] = 0;
        b[1][1] = 1;

        if (printDebug) {
          std::cout << "Before inversion\n";
          nr_explore::printMatrices(a, b);
        }
        nr_explore::gaussj(a, b);
        if (printDebug) {
          std::cout << "After inversion\n";
          nr_explore::printMatrices(a, b);
        }
        
        numType c00 = ( b[0][0] * dfx.X()  + b[1][0] * dfy.X() );
        numType c10 = ( b[0][0] * dfx.Y()  + b[1][0] * dfy.Y() );
        numType c20 = ( b[0][0] * dfx.Z()  + b[1][0] * dfy.Z() );

        numType c01 = ( b[0][1] * dfx.X()  + b[1][1] * dfy.X() );
        numType c11 = ( b[0][1] * dfx.Y()  + b[1][1] * dfy.Y() );
        numType c21 = ( b[0][1] * dfx.Z()  + b[1][1] * dfy.Z() );

        if (printDebug) {
          std::cout << std::fixed << std::setprecision(6) << std::setw(5);

          std::cout << "(JtJ)^{-1} * Jt is \n";
          std::cout << c00 << "\t" << c10 << "\t" << c20 << "\n";
          std::cout << c01 << "\t" << c11 << "\t" << c21 << "\n";

          std::cout << std::defaultfloat;
        }

        numType step_x = c00 * f.X() + c10 * f.Y() + c20 * f.Z(); 
        numType step_y = c01 * f.X() + c11 * f.Y() + c21 * f.Z(); 

        if (printDebug) {
          std::cout << "step_x is " << step_x << "\n";
          std::cout << "step_y is " << step_y << "\n";
        }

        inputT step(step_x, step_y);
        return step;
      }
  };

  Coords2D<double> intersect3DCurves(
    Curve& c1,
    Curve& c2,
    RangeChecker2D<Coords2D<double>>& rangeChecker,
    Coords2D<double>& start
  ) {
    bool printDebug = false;

    using InputType = Coords2D<double>;
    using RangeCheckerType = RangeChecker2D<InputType>;
    using OutputType = Coords3D<double>;
    using FunctionType = CurveDifference3D;
    using StepFinderType = StepFinder2D3D<double, InputType, OutputType, FunctionType>;
    using ConvergenceType = ConvergenceChecker2D<double, InputType>;

    StepFinderType stepFinder;

    ConvergenceType convergence(1e-6);

    // function representing the difference between two lines
    // expect (immediate) convergence to their intersection
    CurveDifference3D diff(c1, c2);

    InputType newtonResultl1l2 = newtonRaphson<
      double,
      InputType,
      RangeCheckerType,
      OutputType,
      FunctionType,
      StepFinderType,
      ConvergenceType
    >(
      diff, 
      start,
      rangeChecker,
      stepFinder,
      convergence, 
      100
    );

    if (printDebug) {
      std::cout << "iteration starting from " <<  start.X() << ", " << start.Y() << " converged to " 
        << newtonResultl1l2.X() << ", " << newtonResultl1l2.Y() << "\n";
    }

    Point onC1 = c1.evaluate(newtonResultl1l2.X());
    Point onC2 = c2.evaluate(newtonResultl1l2.Y());

    if (printDebug) {
      std::cout << "point on c1 at " << newtonResultl1l2.X() << " = " 
        <<  onC1.X() << ", " << onC1.Y() << ", " << onC1.Z() << "\n";
      std::cout << "point on c2 at " << newtonResultl1l2.Y() << " = " 
        <<  onC2.X() << ", " << onC2.Y() << ", " << onC2.Z() << "\n";
    }

    return newtonResultl1l2;
  }

  void plotConvergencePattern() {
    bool printDebug = false;

    Circle c(Point(1.0, 1.0, 5.0), 2);
    Line l(Point(2, 1.0, 5.0), Vector(0, 1, 0));

    int displaySize = 2;

    const int NUM_I = 400;
    const int NUM_J = 400;
    const float LOW_X = - 3 * M_PI;
    const float HIGH_X = 3 * M_PI;
    const float LOW_Y = -10.5;
    const float HIGH_Y = 10.5;

    addCurveToView(c, LOW_X, HIGH_X, RED, 2);
    addCurveToView(l, LOW_Y, HIGH_Y, GREEN, 2);

    std::map<ColorPatch2D<double>, std::vector<Coords2D<double>>> foundSolutions;
    double accuracy_tolerance = 1e-3;


    RangeChecker2D<Coords2D<double>> rangeChecker(
      Coords2D<double>(LOW_X, LOW_Y),
      Coords2D<double>(HIGH_X, HIGH_Y)
    );

    for (int i = 0; i < NUM_I; i++) {
      for (int j = 0; j < NUM_J; j++) {

        Coords2D<double> start(
          LOW_X + (HIGH_X - LOW_X) / NUM_I * i,
          LOW_Y + (HIGH_Y - LOW_Y) / NUM_J * j
        );

        if (printDebug) {
          std::cout << "explore from start(" << start.X() << ", " << start.Y() << ")\n";
        }

        Coords2D<double> newtonResult(
          std::numeric_limits<double>::quiet_NaN(), 
          std::numeric_limits<double>::quiet_NaN()
        );
        try {
          newtonResult = intersect3DCurves(
            c,
            l,
            rangeChecker,
            start
          );
        } catch (const std::runtime_error& e) {
          std::cerr << "intersect3DCurves threw an error: " << e.what() << std::endl;
        } catch (const char* msg) {
          if (std::string(msg).compare("Singular Matrix")) {            
            // don't print an error
            if (printDebug) {
              std::cerr << "Singular Matrix\n";
            }
          } else {
            std::cerr << "intersect3DCurves threw an error: " << msg << std::endl;
          }
        } catch (...) {
          std::cerr << "intersect3DCurves threw an Unknown exception caught!" << std::endl;
        }
        bool addedToMap = false;
        for (auto& kv : foundSolutions) {
          auto& key = kv.first;
          if (patchMatchesCoords(key, newtonResult, accuracy_tolerance * accuracy_tolerance)) {
            if (printDebug) {
              std::cout << "add (" << start.X() << ", " << start.Y() << ") to existing collection for key (" << key.X() << ", " << key.Y() << ")\n";
            }
            kv.second.push_back(start);
            addedToMap = true;
            break;
          }
        }
        if (!addedToMap) {
          // add a new key value pair to the map
          if (printDebug) {
            std::cout << "create new collection for newtonResult (" << newtonResult.X() << ", " << newtonResult.Y() << ")\n";
          }

          if (foundSolutions.size() < colors.size()) {
            ColorPatch2D<double> newPatch(
              newtonResult,
              colors[foundSolutions.size()]
            );
            foundSolutions[
              newPatch
            ] = {start};
          } else {
            std::cout << "Ran out of colors to plot, " 
              << foundSolutions.size() << ">=" 
              << colors.size() << "\n";
          }
        }
      }
    }

    addPtsToView(foundSolutions, displaySize);
  }

  void intersectCurves3D() {

    Circle c(Point(0.0, 0.0, 5.0), 2);
    Line l(Point(1.0, 1.0, 5.0), Vector(0, 1, 0));
    Coords2D<double> start = Coords2D<double>(0.3, 0.1);

    addCurveToView(c, 0.0, 2*M_PI, RED, 2);
    addCurveToView(l, 0.0, 2*M_PI, GREEN, 2);


    RangeChecker2D<Coords2D<double>> rangeChecker(
      Coords2D<double>(-M_PI, -100.0), // lower s, lower t
      Coords2D<double>(M_PI, 100.0)    // upper s, upper t
    );

    Coords2D<double> result = intersect3DCurves(
      c,
      l,
      rangeChecker,
      start
    );

    Point onC1 = c.evaluate(result.X());

    addPointToView(onC1, YELLOW, 7);
  }

  void testNewtonRaphson3Dinput(){
    try {
      Line xaxis(Point(0,0,0), Vector(1,0,0));
      Line yaxis(Point(0,0,0), Vector(0,1,0));
      Line zaxis(Point(0,0,0), Vector(0,0,1));
  
      //addCurveToView(xaxis, 0.0, 1.0, RED, 2);
      //addCurveToView(yaxis, 0.0, 1.0, GREEN, 2);
      //addCurveToView(zaxis, 0.0, 1.0, BLUE, 2);

      //intersectCurves3D();

      plotConvergencePattern();

    } catch (const std::runtime_error& e) {
      std::cerr << "in testNewtonRaphson3Dinput got an error: " << e.what() << std::endl;
    } catch (const char* msg) {
      std::cerr << "in testNewtonRaphson3Dinput got an error: " << msg << std::endl;
    } catch (...) {
      std::cerr << "in testNewtonRaphson3Dinput got an Unknown exception caught!" << std::endl;
    }
  }
}
