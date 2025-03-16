#include <boost/optional/optional.hpp>
#include <iostream>
#include <string>

namespace boost_data_types {

  void optionalExample() {
    boost::optional<int> opt;
    if (opt) {
      std::println("{} Value: {}", __func__, *opt);
    } else {
      std::println("{} No Value", __func__);
    }

    opt = 10;
    if (opt) {
      std::println("{} Value: {}", __func__, *opt);
    } else {
      std::println("{} No Value", __func__);
    }
  }

}

