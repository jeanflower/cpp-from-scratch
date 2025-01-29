import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js'

// This threejs file extracts data from a file output/view_data.json
// It expects to find this json structure
// {
//   "ptsGps": [
//     {
//     "color": "0x00ff00",
//     "pts": [
//       {"x": 0, "y": 0, "z": 0},
//       {"x": 1, "y": 1, "z": 1},
//       ... // add more points here
//     ], // end of points array
//     }, // end of ptsGp object
//     ... // add more ptsGps here
//   ], // end of ptsGps array
//   "linesObjs": [// multiple polylines can be included
//     {
//       "color": "0x00ff00",
//       "vxs": 
//         [ 
//           {"x": 0, "y": 0, "z": 0}, // start of first line segment
//           {"x": 1, "y": 1, "z": 1}, // end of first line segment and start of second
//           ...
//                                   // for closed lines, first and last point should match
//         ],
//     }
//     ... // add more lineObjs here
//   ] // end of lineObjs array
// }
// The data is used to create points and polylines in the scene.


// Set up the scene, camera, and renderer
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
const renderer = new THREE.WebGLRenderer();
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Add orbit controls for interactivity
new OrbitControls(camera, renderer.domElement)

// Parse the XYZ data from the JSON file
fetch("output/view_data.json")
    .then(response => 
      {
        return response.json();
      }
    )
    .then(data => {

        // how to display a given ptsObj (which contains pts data and color)
        const addPointsToScene = (ptsObj) => {
          console.log('processing a ptsObj');
          // Create a geometry to hold points
          const geometry = new THREE.BufferGeometry();
          // make a flattened map of the coordinates
          const points = new Float32Array(ptsObj.pts.flatMap(d => [d.x, d.y, d.z]));
          // console.log(`point data is ${points}`);
          geometry.setAttribute("position", new THREE.BufferAttribute(points, 3));
          // Create a material for the points
          const material = new THREE.PointsMaterial({ color: parseInt(ptsObj.color, 16), size: 0.05 });
          const pointCloud = new THREE.Points(geometry, material);
          // Add the points to the scene
          scene.add(pointCloud);
        }

        // for each ptsOj, display it
        data.ptsGps.map(addPointsToScene);

        // how to display a given linesObj (which contains vxs data and color)
        const addPolylineToScene = (linesObj) => { // TODO refactor to share code with addPointsToScene
          // Create a geometry to hold line segment vertices
          const geometry = new THREE.BufferGeometry();
          // Create a flattened array of vertex coordinates for the line segments
          const vertices = new Float32Array(linesObj.vxs.flatMap(d => [d.x, d.y, d.z]));
          // console.log(`polyline data is ${vertices}`);
          geometry.setAttribute("position", new THREE.BufferAttribute(vertices, 3));
      
          // Create a material for the line segments
          const material = new THREE.LineBasicMaterial({ color: parseInt(linesObj.color, 16), linewidth: 3 });
      
          // Create the line segments and add them to the scene
          const lineSegments = new THREE.Line(geometry, material);
          scene.add(lineSegments);
        };
        // for each linesObj, display it
        data.linesObjs.map(addPolylineToScene);

        // Animate the scene
        animate();
    });

// Set the camera position
camera.position.z = 20; // TODO set camera for different data

// Animation loop
function animate() {
    requestAnimationFrame(animate);
    // controls.update();
    renderer.render(scene, camera);
}
