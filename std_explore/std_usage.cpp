#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include <list>

namespace std_data_types {

  void ptr_example() {

    // for testing leak detection
    //const int* leaking_ptr = new int(10);

    const std::unique_ptr<int> sp1(new int(10));
    std::cout << __func__ << " " << "Value: " << *sp1 << std::endl;
    // when the unique_ptr goes out of scope, the memory is automatically deallocated

    int* co_owned_ptr = new int(20);
    const std::shared_ptr<int> sp2(co_owned_ptr);
    std::cout << __func__ << " " << "sp2 Reference count: " << sp2.use_count() << std::endl;
    {
      const std::shared_ptr<int> sp3(sp2);
      std::cout << __func__ << " " << "Value: " << *sp3 << std::endl;
      std::cout << __func__ << " " << "sp2 Reference count: " << sp2.use_count() << std::endl;
    }
    std::cout << __func__ << " " << "sp2 Reference count: " << sp2.use_count() << std::endl;
    // when the last shared_ptr goes out of scope, the memory is automatically deallocated
  }


void string_example() {
    std::string str = "Standard C++ Libraries";

    // Convert to lowercase
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);

    std::cout << __func__ << " " << "Lowercase string: " << str << std::endl;

    // Check if the string contains "boost"
    if (str.find("standard") != std::string::npos) {
        std::cout << __func__ << " " << "String contains 'standard'!" << std::endl;
    }
  }

  void collections_example() {
    /*
    When to Use Which Collection:

    std::vector: 
    Use when you need random access to elements and the ability to dynamically resize the collection.

    std::list: 
    Use when you need efficient insertions/deletions in the middle of the collection but don't need fast random access.

    std::deque: 
    Use when you need efficient insertions and deletions at both ends of the collection.

    std::set and std::map: 
    Use when you need to store unique elements (or key-value pairs) in sorted order.

    std::unordered_set and std::unordered_map: 
    Use when you need efficient lookups, insertions, and deletions without caring about sorting 
    order.

    std::queue and std::stack: 
    Use when you need to maintain order in the collection and need efficient access to either 
    the front (queue) or the top (stack).

    std::array: 
    Use when the collection size is known at compile time and you want the performance benefits of a 
    fixed-size array but still want some flexibility over a raw array.
    */
    const int BIG_NUMBER = 800000;

    auto fillInsert = [BIG_NUMBER](auto& coll) {
      for (int i = 0; i < BIG_NUMBER; i++) {
        auto it = coll.begin();
        if (i > 2) {
          it = std::next(coll.begin(), 2);
        }
        coll.insert(it, i);
      }
    };

    /////////////////////////
    /*   Vector utilities  */
    /////////////////////////
    auto makeNewVector = []() {
      return std::vector<int>();
    };

    auto fillAndQueryVector = [BIG_NUMBER](auto fillFunc, auto queryFunc) {
      std::vector<int> vec = {};
      fillFunc(vec);
      return queryFunc(vec);
    };

    // N.B. BIG_NUMBER is provided in a capture clause
    auto fillVecPush = [BIG_NUMBER](std::vector<int>& vec) {
      for (int i = 0; i < BIG_NUMBER; i++) {
        vec.push_back(i);
      }
    };

    auto queryVecSum = [BIG_NUMBER](std::vector<int>& vec) {
      std::uint64_t sum = 0;
      for (int i = 0; i < BIG_NUMBER; i++) {
        sum += vec[i];
      }
      return sum;
    };

    auto queryVecLength = [BIG_NUMBER](std::vector<int>& vec) {
      return vec.size();
    };


    /////////////////////////
    /*   List utilities    */
    /////////////////////////
    auto makeNewList = []() {
      return std::list<int>();
    };

    auto fillAndQueryList = [BIG_NUMBER](auto fillFunc, auto queryFunc) {
      std::list<int> lst = {};
      fillFunc(lst);
      return queryFunc(lst);
    };

    auto fillListPush = [BIG_NUMBER](std::list<int>& lst) {
      for (int i = 0; i < BIG_NUMBER; i++) {
        lst.push_back(i);
      }
    };

    auto queryListSum = [BIG_NUMBER](std::list<int>& lst) {
      std::uint64_t sum = 0;
      for (const int elem : lst) {
        sum += elem;
      }
      return sum;
    };

    auto queryListLength = [BIG_NUMBER](std::list<int>& lst) {
      return lst.size();
    };

    /////////////////////////
    /*   Deque utilities   */
    /////////////////////////
    auto makeNewDeque = []() {
      return std::deque<int>();
    };

    auto fillAndQueryDeque = [BIG_NUMBER](auto fillFunc, auto queryFunc) {
      std::deque<int> dq = {};
      fillFunc(dq);
      return queryFunc(dq);
    };

    auto fillDequePush = [BIG_NUMBER](std::deque<int>& dq) {
      for (int i = 0; i < BIG_NUMBER; i++) {
        dq.push_back(i);
      }
    };

    auto queryDequeSum = [BIG_NUMBER](std::deque<int>& dq) {
      std::uint64_t sum = 0;
      for (const int elem : dq) {
        sum += elem;
      }
      return sum;
    };

    auto queryDequeLength = [BIG_NUMBER](std::deque<int>& dq) {
      return dq.size();
    };

    auto timeWork = [](
      auto makeNewCollection,
      auto fillCollection, 
      auto queryCollection,
      const std::string& message
    ) {
      auto timeOnce = [makeNewCollection, fillCollection, queryCollection](
      ) { 
        // Start timer
        auto start = std::chrono::high_resolution_clock::now();

        auto coll = makeNewCollection();
        fillCollection(coll);
        const std::uint64_t result = queryCollection(coll);

        // End timer
        auto end = std::chrono::high_resolution_clock::now();
        // Calculate duration
        std::chrono::duration<double> duration = end - start;
        return duration.count();
      };
      auto d1 = timeOnce();
      auto result = d1;

      std::cout << __func__ << " " << "Time taken is " << result << "s using " << message << std::endl;
    };

    std::cout << __func__ << "------------ push-sum -- deque has faster push" << std::endl;
    timeWork(makeNewVector, fillVecPush, queryVecSum, "vector-push-sum");
    timeWork(makeNewList, fillListPush, queryListSum, "list-push-sum");
    timeWork(makeNewDeque, fillDequePush, queryDequeSum, "deque-push-sum");

    std::cout << __func__ << "------------ push-length -- deque has faster push" << std::endl;
    timeWork(makeNewVector, fillVecPush, queryVecLength, "vector-push-length");
    timeWork(makeNewList, fillListPush, queryListLength, "list-push-length");
    timeWork(makeNewDeque, fillDequePush, queryDequeLength, "deque-push-length");

    std::cout << __func__ << "------------ insert-sum -- list has faster insert - vector is terrible!" << std::endl;
    timeWork(makeNewVector, fillInsert, queryVecSum, "vector-insert-sum");
    timeWork(makeNewList, fillInsert, queryListSum, "list-insert-sum");
    timeWork(makeNewDeque, fillInsert, queryDequeSum, "deque-insert-sum");

    std::cout << __func__ << "------------ insert-length -- list has faster insert - vector is terrible!" << std::endl;
    timeWork(makeNewVector, fillInsert, queryVecLength, "vector-insert-length");
    timeWork(makeNewList, fillInsert, queryListLength, "list-insert-length");
    timeWork(makeNewDeque, fillInsert, queryDequeLength, "deque-insert-length");
  }
}
