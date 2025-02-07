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
    std::optional<std::reference_wrapper<Vector>> first_derivative, 
    std::optional<std::reference_wrapper<Vector>> second_derivative
  ) const {
    Point p{radius * cos(t), radius * sin(t)};      

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
    Circle c(5.0);
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
    const Circle c(12.0);

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
}
