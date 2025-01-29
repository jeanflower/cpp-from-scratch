import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js'
import { MeshLine, MeshLineMaterial } from "three.meshline";

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
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Add orbit controls for interactivity
new OrbitControls(camera, renderer.domElement)

function fitCameraToScene() {
  const box = new THREE.Box3();
  const tempBox = new THREE.Box3();
  const center = new THREE.Vector3();

  // Expand bounding box to include all objects
  scene.traverse((object) => {
    // console.log(`object is ${JSON.stringify(object)}`);
    if (object.isMesh || object.isObject3D) {
      object.updateWorldMatrix(true, true); // Ensure world transformation is applied

      // Ensure the object has geometry (some Object3Ds may not)
      if (object.geometry) {
        tempBox.setFromObject(object);
        box.expandByPoint(tempBox.min);
        box.expandByPoint(tempBox.max);
      }
    }
  });

  // Check if the scene contains any objects
  if (box.isEmpty()) {
      console.warn("Scene is empty or contains no meshes.");
      return;
  }

  box.getCenter(center);
  const size = box.getSize(new THREE.Vector3());

  // Determine the camera distance required
  const maxDim = Math.max(size.x, size.y, size.z);
  const fov = THREE.MathUtils.degToRad(camera.fov);
  const cameraDistance = maxDim / (2 * Math.tan(fov / 2));

  if (camera.isPerspectiveCamera) {
      camera.position.set(center.x, center.y, center.z + cameraDistance * 1.5);
  } else if (camera.isOrthographicCamera) {
      camera.left = -size.x / 2;
      camera.right = size.x / 2;
      camera.top = size.y / 2;
      camera.bottom = -size.y / 2;
      camera.near = -size.z;
      camera.far = size.z * 2;
      camera.updateProjectionMatrix();
  }

  camera.lookAt(center);
}

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
          const col = ptsObj.color;
          // console.log(`col = ${col}`);

          // Create a geometry to hold points
          const geometry = new THREE.BufferGeometry();
          // make a flattened map of the coordinates
          const points = new Float32Array(ptsObj.pts.flatMap(d => [d.x, d.y, d.z]));
          // console.log(`point data is ${points}`);
          geometry.setAttribute("position", new THREE.BufferAttribute(points, 3));
          // Create a material for the points
          const material = new THREE.PointsMaterial({ 
            color: col, 
            size: 0.1
          });
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

          const line = new MeshLine();
          line.setPoints(vertices); // Example points

          const material = new MeshLineMaterial({
            color: new THREE.Color(linesObj.color),
            lineWidth: 0.05,  // Works properly
          });

          const mesh = new THREE.Mesh(line, material);
          scene.add(mesh);
        };
        // for each linesObj, display it
        data.linesObjs.map(addPolylineToScene);

        fitCameraToScene();

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
