#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "boost_explore/boost_usage.hpp"
#include "std_explore/std_usage.hpp"
#include "geom_explore/geomAPI_usage.hpp"
#include "geom_explore/geom.hpp"
#include <iostream>
#include <string>

TEST_CASE("Sample Test", "[sample]") {
    REQUIRE(1 + 1 == 2);  // This test will pass
    //REQUIRE(1 + 1 == 3);  // This test will fail
}
TEST_CASE("Boost", "[coverage]") {
    boost_data_types::optionalExample();
}
TEST_CASE("std", "[coverage]") {
    std_data_types::stringExample();
    std_data_types::ptrExample();
    std_data_types::collectionsExample();
}
TEST_CASE("cascade", "[coverage]") {
    geomAPI_examples::sphereExample();
    geomAPI_examples::torusExample();
}
