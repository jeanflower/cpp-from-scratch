#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "boost_explore/boost_usage.hpp"
#include "std_explore/std_usage.hpp"
#include "geom_explore/geomAPI_usage.hpp"
#include "geom_explore/geom.hpp"
#include <iostream>
#include <string>
#include "nr_explore/gaussj.hpp"

TEST_CASE("Sample Test", "[sample]") {
    REQUIRE(1 + 1 == 2);  // This test will pass

    boost_data_types::optionalExample();

    std_data_types::stringExample();
    std_data_types::ptrExample();
    std_data_types::collectionsExample();
    
    geomAPI_examples::sphereExample();
    geomAPI_examples::torusExample();

    geom_examples::nurbsExample();
    const std::string message = geom_examples::nurbsPerformanceExample();
    std::cout << message;

    geom_examples::writeGeometryToJSON();

    nr_explore::testGaussj();

    geom_examples::circleExample();

}
