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

  Line::Line(
    Point start,
    Vector direction
  ) : start(start), direction(direction) {}

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

  Circle::Circle(
    Point origin,
    double r
  ) : origin(origin), radius(r) {}

  Point Circle::evaluate(
    double t, 
    std::optional<std::reference_wrapper<Vector>> first_derivative, 
    std::optional<std::reference_wrapper<Vector>> second_derivative
  ) const {
    Point p{
      origin.X() + radius * cos(t), 
      origin.Y() + radius * sin(t),
      origin.Z()
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
  
  void circleExample() {
    Circle c(Point(0.0, 0.0, 0.0), 5.0);
    double t = 1.0;

    Vector first, second;
    Point p = c.evaluate(t, first, second);

    std::cout << "Point: (" << p.x << ", " << p.y << ")\n";
    std::cout << "First Derivative: (" << first.x << ", " << first.y << ")\n";
    std::cout << "Second Derivative: (" << second.x << ", " << second.y << ")\n";
  }

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

  // can be used for testing derivatives at any order
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

  template <class T>
  struct ComplexNumber {
    T x, y;

    // Constructor
    ComplexNumber(T x = 0, T y = 0);

    // Getter functions
    T X() const;
    T Y() const;

    T len_sq() const;
  };

  template <class T>
  ComplexNumber<T>::ComplexNumber(T x, T y) : x(x), y(y) {}

  template <class T>
  T ComplexNumber<T>::X() const { return x; }
  template <class T>
  T ComplexNumber<T>::Y() const { return y; }

  template <class T>
  T ComplexNumber<T>::len_sq() const { return x * x + y * y; }

  template <class T>
  ComplexNumber<T> operator+(const ComplexNumber<T>& a, const ComplexNumber<T>& b) {
    return ComplexNumber(a.x + b.x, a.y + b.y);
  }

  template <class T>
  ComplexNumber<T> operator-(const ComplexNumber<T>& a, const ComplexNumber<T>& b) {
    return ComplexNumber(a.x - b.x, a.y - b.y);
  }

  template <class T>
  ComplexNumber<T> operator*(const ComplexNumber<T>& a, const ComplexNumber<T>& b) {
    return ComplexNumber(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
  }

  template <class T>
  class ComplexFunctionExample {
    public:
      explicit ComplexFunctionExample() {};

      ComplexNumber<T> operator() (const ComplexNumber<T> z) {
        return ComplexNumber<T>(
          z.X() * z.X() * z.X() - 3 * z.X() * z.Y() * z.Y() - 1,
          3 * z.X() * z.X() * z.Y() - z.Y() * z.Y() * z.Y()
        );
      }
      ComplexNumber<T> dfx(const ComplexNumber<T> x) {
        return ComplexNumber<T>(
          3 * x.X() * x.X() - 3 * x.Y() * x.Y(),
          -6 * x.X() * x.Y()
        );
      }
      ComplexNumber<T> dfy(const ComplexNumber<T> x) {
        return ComplexNumber<T>(
          6 * x.X() * x.Y(),
          3 * x.X() * x.X() - 3 * x.Y() * x.Y()
        );
      }
  };

  template <class T, class F>
  ComplexNumber<T> complex_newt(
    F &funcd, 
    ComplexNumber<T>& guess,
    const T lower_bound_x, 
    const T upper_bound_x, 
    const T lower_bound_y, 
    const T upper_bound_y, 
    const T accuracy_tolerance,
    const int MAX_ITERATIONS = 200
  ) {
    const double accuracy_tolerance_sq = accuracy_tolerance * accuracy_tolerance;
    ComplexNumber<T> rtn = guess;
    for (int j = 0; j < MAX_ITERATIONS; j++) {
      ComplexNumber<T> f = funcd(rtn);
      ComplexNumber<T> dfx = funcd.dfx(rtn);
      ComplexNumber<T> dfy = funcd.dfy(rtn);

      double det = dfx.X() * dfy.Y() - dfx.Y() * dfy.X();
      if(det == 0) {
          throw std::runtime_error("Zero derivative determinant in complex_newt");
      }
      
      double step_x = ( f.X() * dfy.Y() - f.Y() * dfx.Y() ) / det;
      double step_y = ( f.Y() * dfx.X() - f.X() * dfy.X() ) / det;
      ComplexNumber<T> step(step_x, step_y);

      rtn = rtn - step;

      // std::cout << "j is " << j << " and rtn is " << rtn.X() << ", " << rtn.Y() << "\n";
      if ((lower_bound_x - rtn.X()) * (rtn.X() - upper_bound_x) < 0.0 ||
        (lower_bound_y - rtn.Y()) * (rtn.Y() - upper_bound_y) < 0.0)
        throw("Jumped out of range in complex_newt");
      if (step.len_sq() < accuracy_tolerance_sq) {
        // std::cout << "Converged to " << rtn.X() << ", " << rtn.Y() << " after " << j << " iterations\n";
        return rtn; 
      }
    }
    throw("Maximum number of iterations exceeded in complex_newt");
  }

  template <class T>
  bool complexNRConverges(
    ComplexFunctionExample<T>& f,
    ComplexNumber<T>& guess
  ) {
    try {
      ComplexNumber<T> start = guess;

      //std::cout << "newtonResult on f from guess " << start.X() << " + i * " << start.Y() << "\n";
      ComplexNumber<T> newtonResult = complex_newt<T>(
        f, 
        start,
        -100.0, 100.0, 
        -100.0, 100.0, 
        1e-10,
        100
      );
      //std::cout << newtonResult.X() << " + i * " << newtonResult.Y() << "\n";

      if (abs(newtonResult.X() - 1.0) < 0.1 && abs(newtonResult.Y() - 0.0) < 0.1) {
        // we converged on the root of interest
        return true;
      } else {
        return false;
      }

    } catch (const char* msg) {
      //std::cerr << "Error: " << msg << std::endl;
    } catch (...) {
      //std::cerr << "Unknown exception caught!" << std::endl;
    }
    return false;
  }

  void fractal() {

    std::vector<Point> ptsConverged;
    std::vector<Point> ptsNotConverged;

    ComplexFunctionExample f = ComplexFunctionExample<long double>();
    const int NUM_I = 100;
    const int NUM_J = 100;
    const long double LOW_X = -2.0;
    const long double HIGH_X = 2.0;
    const long double LOW_Y = -2.0;
    const long double HIGH_Y = 2.0;
    for (int i = 0; i < NUM_I; i++) {
      for (int j = 0; j < NUM_J; j++) {

        ComplexNumber<long double> start(
          LOW_X + (HIGH_X - LOW_X) / NUM_I * i,
          LOW_Y + (HIGH_Y - LOW_Y) / NUM_J * j 
        );

        bool converged = complexNRConverges(f, start);
        //std::cout << start.X() << " + i * " << start.Y();
        if (converged) {
          //std::cout << " Converged\n";
          ptsConverged.push_back(Point(start.X(), start.Y(), 0));
          // std::cout << "+";
        } else {
          //std::cout << " Did not converge\n";
          ptsNotConverged.push_back(Point(start.X(), start.Y(), 0));
          // std::cout << "-";
        }
      }
    }

    std::vector<PtCollection> ptsColls;
    PtCollection ptsObjConverged ={
      .displaySize = 1,
      .color = BLACK,
      .pts = ptsConverged,
      .isLine = false
    };
    ptsColls.push_back(ptsObjConverged);
    PtCollection ptsObjNotConverged ={
      .displaySize = 1,
      .color = WHITE,
      .pts = ptsNotConverged,
      .isLine = false
    };
    ptsColls.push_back(ptsObjNotConverged);
    // std::cout << "adding " << ptsObj.pts.size() << " to view\n";
    addGeometryToView(ptsColls);
  }

}
