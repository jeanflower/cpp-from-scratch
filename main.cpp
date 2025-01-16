#include <boost/algorithm/string.hpp>
#include <iostream>
#include <string>

using namespace std;

int main() {
    std::string str = "Boost C++ Libraries";
    
    // Convert string to lowercase using Boost
    boost::algorithm::to_lower(str);

    std::cout << "Lowercase string: " << str << std::endl;

    // Check if the string contains "boost"
    if (boost::algorithm::contains(str, "boost")) {
        std::cout << "String contains 'boost'!" << std::endl;
    }

    return 0;
}
