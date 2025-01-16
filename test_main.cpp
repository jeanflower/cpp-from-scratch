#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>

TEST_CASE("Sample Test", "[sample]") {
    REQUIRE(1 + 1 == 2);  // This test will pass
}
