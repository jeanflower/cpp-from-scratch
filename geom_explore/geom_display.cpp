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

namespace geom_examples {

  std::vector<geom_examples::PtCollection> viewables;

  void addGeometryToView(const std::vector<PtCollection>& ptColls) {
    for (auto x : ptColls) {
      viewables.push_back(x);
    }
  }

  void writeGeometryToJSON(
  ) {
    try {
      // Write this data to a json file for viewing to pick up
      std::ofstream debugFile("viewer/public/output/view_data.json");
      if (!debugFile) {
          std::cerr << "Error opening file for writing.\n";
          return;
      }

      std::vector<PtCollection> nonLineColls;
      std::vector<PtCollection> isLineColls;
      std::copy_if(
        viewables.begin(), 
        viewables.end(), 
        std::back_inserter(nonLineColls), 
        [](const PtCollection& elt) { return !elt.isLine; }
      );
      std::copy_if(
        viewables.begin(), 
        viewables.end(), 
        std::back_inserter(isLineColls), 
        [](const PtCollection& elt) { return elt.isLine; }
      );

      // we'll redirect cout to write to json file
      std::streambuf* coutBuf = std::cout.rdbuf(); // Save the original buffer

      std::cout.rdbuf(debugFile.rdbuf());  // Redirect std::cout to the file
      std::cout << "{";
      std::cout << "\"ptsGps\":[\n";
      for (size_t i = 0; i < nonLineColls.size(); ++i) {
        auto pts = nonLineColls[i].pts;
        auto col = nonLineColls[i].color;
        auto displaySize = nonLineColls[i].displaySize;
        std::cout << "{";
        std::cout << "\"color\":" << col << ",\n";
        std::cout << "\"displaySize\":" << displaySize << ",\n";
        std::cout << "\"pts\":[\n";
        for (size_t j = 0; j < pts.size(); ++j) {
          const Point& pt = pts[j];
          std::cout 
            << "  { \"x\": " << pt.X() << ", \"y\": " << pt.Y() << ", \"z\": " << pt.Z() << "}";
          if (j != pts.size() - 1) {
            std::cout << ",";
          }
          std::cout << "\n"; // new point on a new line of the json file
        }
        std::cout << "]\n"; // end of array of points
        std::cout << "}\n";  // end of pts object
        if (i != nonLineColls.size() - 1) {
          std::cout << ",";
        }
      }
      std::cout << "],\n"; // end of array of pts objects
      std::cout << "\"linesObjs\":[\n";

      for (size_t polyLineIndex = 0; polyLineIndex < isLineColls.size(); ++polyLineIndex) {
        auto lineColl = isLineColls[polyLineIndex];
        auto col = lineColl.color;
        auto vxs = lineColl.pts;
        auto displaySize = lineColl.displaySize;
        std::cout << "{";
        std::cout << "\"color\":" << col << ",\n";
        std::cout << "\"displaySize\":" << displaySize << ",\n";
        std::cout << "\"vxs\":";
        std::cout << "[\n";
        for (size_t segIndex = 0; segIndex < vxs.size(); ++segIndex) {
          const Point& pt = vxs[segIndex];
          std::cout 
            << "  { \"x\": " << pt.X() << ", \"y\": " << pt.Y() << ", \"z\": " << pt.Z() << "}";
          if (segIndex != vxs.size() - 1) {
            std::cout << ",";
          }
          std::cout << "\n"; // new point on a new line of the json file
        }
        std::cout << "]\n"; // end of vxs array
        std::cout << "}\n"; // end of polyLine object

        if (polyLineIndex != isLineColls.size() - 1) {
          std::cout << ",";
        }
        std::cout << "\n";
      }
      std::cout << "]}";
      // Restore std::cout to its original state
      std::cout.rdbuf(coutBuf);
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Unknown exception caught!" << std::endl;
    }
  }

  void addCurveToView(
    const Curve& c,
    double start_t,
    double end_t,
    uint32_t color,
    int displaySize
  ) {
    // evaluate at a set of points and display
    std::vector<Point> ptsToDisplay;
    const int NUM_SAMPLES = 3000;
    for (int i = 0; i < NUM_SAMPLES; i++) {
      const double t = start_t + i * (end_t - start_t) / (NUM_SAMPLES - 1);
      const Point& result = c.evaluate(t);
      ptsToDisplay.push_back(Point(result.x, result.y, result.z));
      // std::cout << t << ", " 
      // << result.x << ", " << result.y << ", " << result.z << std::endl;
    }
    std::vector<PtCollection> ptsColls;
    PtCollection ptsObj ={
      .displaySize = displaySize,
      .color = color,
      .pts = ptsToDisplay,
      .isLine = false
    };
    ptsColls.push_back(ptsObj);
    std::cout << "adding " << ptsObj.pts.size() << " to view\n";
    addGeometryToView(ptsColls);
  }

  void addPointToView(
    const Point& p,
    uint32_t color,
    int displaySize 
  ) {
    // evaluate at a set of points and display
    std::vector<Point> ptsToDisplay;
    ptsToDisplay.push_back(p);
    std::vector<PtCollection> ptsColls;
    PtCollection ptsObj ={
      .displaySize = displaySize,
      .color = color,
      .pts = ptsToDisplay,
      .isLine = false
    };
    ptsColls.push_back(ptsObj);
    std::cout << "adding ()" << p.X() << ", " << p.Y() << ", " << p.Z() << ") to view\n";
    addGeometryToView(ptsColls);
  }
}
