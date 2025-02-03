#include "boost_explore/boost_usage.hpp"
#include "std_explore/std_usage.hpp"
#include "geom_explore/geomAPI_usage.hpp"
#include "geom_explore/geom.hpp"
#include <iostream>
#include <string>
#include "nr_explore/gaussj.hpp"

int main() {
    //boost_data_types::optional_example();

    //std_data_types::string_example();
    //std_data_types::ptr_example();
    //std_data_types::collections_example();
    
    //geomAPI_examples::sphere_example();
    //geomAPI_examples::torus_example();

    //geom_examples::nurbs_example();
    //const std::string message = geom_examples::nurbs_performance_example();
    //std::cout << message;

    //geom_examples::writeGeometryToJSON();

    nr_explore::testGaussj();

    return 0;
}
