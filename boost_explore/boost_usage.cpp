#include <boost/optional/optional.hpp>
#include <iostream>
#include <string>

namespace boost_data_types {

  void optionalExample() {
    boost::optional<int> opt;
    if (opt) {
      std::cout << __func__ << " " << "Value: " << *opt << std::endl;
    } else {
      std::cout << __func__ << " " << "No value" << std::endl;
    }

    opt = 10;
    if (opt) {
      std::cout << __func__ << " " << "Value: " << *opt << std::endl;
    } else {
      std::cout << __func__ << " " << "No value" << std::endl;
    }
  }

}

