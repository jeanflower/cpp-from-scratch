#include <google/cloud/functions/framework.h>
#include <cstdlib>
#include "geom_explore/geom.hpp"

namespace gcf = ::google::cloud::functions;

auto hello_world_http() {
  return gcf::MakeFunction([](gcf::HttpRequest const& /*request*/) {
    std::string greeting = "Hello ";
    auto const* target = std::getenv("TARGET");
    greeting += target == nullptr ? "World" : target;
    greeting += "\n";

    std::string perf_data = geom_examples::nurbs_performance_example();
    greeting += perf_data;
    
    return gcf::HttpResponse{}
        .set_header("Content-Type", "text/plain")
        .set_payload(greeting);
  });
}

int main(int argc, char* argv[]) {
  return gcf::Run(argc, argv, hello_world_http());
}