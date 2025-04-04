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

#include "../rapidjson/document.h"
#include "../rapidjson/writer.h"
#include "../rapidjson/stringbuffer.h"

namespace geom_examples {

  std::vector<geom_examples::PtCollection> viewables;

  void clearView() {
    viewables.clear();
  }
  void addGeometryToView(const std::vector<PtCollection>& ptColls) {
    for (auto x : ptColls) {
      viewables.push_back(x);
    }
  }

  void writeGeometryToJSON(
  ) {
    std::cout << "writing geometry to json\n";
    try {
      // Write this data to a json file for viewing to pick up

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

      std::ofstream debugFile("viewer/public/output/view_data.json");
      if (!debugFile) {
          std::cerr << "Error opening file for writing.\n";
          return;
      }

      rapidjson::Document doc;
      doc.SetObject();
      rapidjson::Document::AllocatorType& allocator = doc.GetAllocator();

      rapidjson::Value ptsGps(rapidjson::kArrayType);
      for (const PtCollection& item : nonLineColls) {
          // Convert item to a Value (depends on its type)
          rapidjson::Value obj(rapidjson::kObjectType); // Create an empty JSON object
    
          obj.AddMember("color", item.color, allocator);
          obj.AddMember("displaySize", item.displaySize, allocator);

          // Convert std::vector<Point> -> JSON Array
          rapidjson::Value ptsArray(rapidjson::kArrayType);
          for (const auto& p : item.pts) {
              rapidjson::Value pointObj(rapidjson::kObjectType);
              pointObj.AddMember("x", p.x, allocator);
              pointObj.AddMember("y", p.y, allocator);
              pointObj.AddMember("z", p.z, allocator);
              ptsArray.PushBack(pointObj, allocator);
          }
          obj.AddMember("pts", ptsArray, allocator); // Add the array to the object
          ptsGps.PushBack(obj, allocator);
      }

      rapidjson::Value linesObjs(rapidjson::kArrayType);
      for (const PtCollection& item : isLineColls) {
          // Convert item to a Value (depends on its type)
          rapidjson::Value obj(rapidjson::kObjectType); // Create an empty JSON object
    
          obj.AddMember("color", item.color, allocator);
          obj.AddMember("displaySize", item.displaySize, allocator);

          // Convert std::vector<Point> -> JSON Array
          rapidjson::Value ptsArray(rapidjson::kArrayType);
          for (const auto& p : item.pts) {
              rapidjson::Value pointObj(rapidjson::kObjectType);
              pointObj.AddMember("x", p.x, allocator);
              pointObj.AddMember("y", p.y, allocator);
              pointObj.AddMember("z", p.z, allocator);

              ptsArray.PushBack(pointObj, allocator);
          }
          obj.AddMember("pts", ptsArray, allocator); // Add the array to the object
          linesObjs.PushBack(obj, allocator);
      }

      rapidjson::Value jsonObj(rapidjson::kObjectType);
      jsonObj.AddMember("ptsGps", ptsGps, allocator);
      jsonObj.AddMember("linesObjs", linesObjs, allocator);

      rapidjson::StringBuffer buffer;
      rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
      jsonObj.Accept(writer); // Accept the writer to generate the JSON
      debugFile << buffer.GetString(); // Write the JSON string to the file
      
      debugFile.flush();
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    } catch (...) {
      std::cerr << "Unknown exception caught! in writeGeometryToJSON" << std::endl;
    }
    std::cout << "finished writing geometry to json\n";
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
    // std::cout << "adding " << ptsObj.pts.size() << " to view\n";
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
    // std::cout << "adding (" << p.X() << ", " << p.Y() << ", " << p.Z() << ") to view\n";
    addGeometryToView(ptsColls);
  }


  // Convert HSL to RGB and return as packed uint32_t (0xRRGGBB)
  uint32_t HSLtoRGB(float h, float s, float l) {
    auto hue2rgb = [](float p, float q, float t) {
      if (t < 0.0f) t += 1.0f;
      if (t > 1.0f) t -= 1.0f;
      if (t < 1.0f / 6.0f) return p + (q - p) * 6.0f * t;
      if (t < 1.0f / 2.0f) return q;
      if (t < 2.0f / 3.0f) return p + (q - p) * (2.0f / 3.0f - t) * 6.0f;
      return p;
    };

    float q = (l < 0.5f) ? l * (1.0f + s) : (l + s - l * s);
    float p = 2.0f * l - q;

    float r = hue2rgb(p, q, h + 1.0f / 3.0f);
    float g = hue2rgb(p, q, h);
    float b = hue2rgb(p, q, h - 1.0f / 3.0f);

    return ((uint8_t)(r * 255) << 16) | ((uint8_t)(g * 255) << 8) | (uint8_t)(b * 255);
  }

  // Generate N muted colors
  std::vector<uint32_t> generateMutedColors(int count) {
    std::vector<uint32_t> colors;
    float saturation = 0.45f;  // Muted but not gray
    float lightness = 0.55f;   // Mid-lightness for good contrast

    for (int i = 0; i < count; ++i) {
      float hue = (i * 1.61803398875f);  // Golden ratio for even distribution
      hue = std::fmod(hue, 1.0f);        // Keep in range [0,1]
      colors.push_back(HSLtoRGB(hue, saturation, lightness));
    }

    return colors;
  }
}
