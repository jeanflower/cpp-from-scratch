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

const aspect = window.innerWidth / window.innerHeight;
const frustumSize = 10; // Controls how much of the scene is visible

const camera = new THREE.OrthographicCamera(
  -frustumSize * aspect,  // Left
   frustumSize * aspect,  // Right
   frustumSize,           // Top
  -frustumSize,           // Bottom
   0.1,                   // Near
   1000                   // Far
);

const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setSize(window.innerWidth, window.innerHeight);
document.body.appendChild(renderer.domElement);

// Add orbit controls for interactivity
// Create and store OrbitControls in a variable
const controls = new OrbitControls(camera, renderer.domElement);
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
  const aspect = window.innerWidth / window.innerHeight;

  if (camera.isPerspectiveCamera) {
    const fov = THREE.MathUtils.degToRad(camera.fov);
    const cameraDistance = maxDim / (2 * Math.tan(fov / 2));
    camera.position.set(center.x, center.y, center.z + cameraDistance * 1.5);
  } 
  else if (camera.isOrthographicCamera) {
    const frustumSize = maxDim * 1.1; // Adjust scale factor to fit scene well
    
    if (aspect >= 1) {
      // Wide screen: X-axis is wider
      camera.left = -frustumSize * aspect / 2;
      camera.right = frustumSize * aspect / 2;
      camera.top = frustumSize / 2;
      camera.bottom = -frustumSize / 2;
    } else {
      // Tall screen: Y-axis is taller
      camera.left = -frustumSize / 2;
      camera.right = frustumSize / 2;
      camera.top = (frustumSize / aspect) / 2;
      camera.bottom = -(frustumSize / aspect) / 2;
    }

    // TODO understand why I needed a +- 50 here
    camera.near = Math.min(box.min.x, box.min.y, box.min.z) - 50;
    camera.far  = Math.max(box.max.x, box.max.y, box.max.z) + 50;

    // console.log(`box = ${JSON.stringify(box)}`);
    // console.log(`center = ${JSON.stringify(center)}`);
    // console.log(`size = ${JSON.stringify(size)}`);
    // console.log(`maxDim = ${JSON.stringify(maxDim)}`);
    // console.log(`camera.left = ${JSON.stringify(camera.left)}`);
    // console.log(`camera.right = ${JSON.stringify(camera.right)}`);
    // console.log(`camera.top = ${JSON.stringify(camera.top)}`);
    // console.log(`camera.bottom = ${JSON.stringify(camera.bottom)}`);
    // console.log(`camera.near = ${JSON.stringify(camera.near)}`);
    // console.log(`camera.far = ${JSON.stringify(camera.far)}`);

    camera.updateProjectionMatrix();
  }

  camera.lookAt(center);

  // **Update OrbitControls target so it doesn't jump when interacting**
  controls.target.copy(center);
  controls.update();
}

// Parse the XYZ data from the JSON file
fetch("output/view_data.json")
    .then(response => 
      {
        return response.json();
      }
    )
    .then(data => {

        // display a ptsObj (which contains pts data and color)
        data.ptsGps.map((ptsObj) => {
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
            size: ptsObj.displaySize,
            sizeAttenuation: false  // displaySize is interpreted as pixels.
          });
          const pointCloud = new THREE.Points(geometry, material);
          // Add the points to the scene
          scene.add(pointCloud);
        });

        // display a linesObj (which contains vxs data and color)
        data.linesObjs.map((linesObj) => {
          // Create a geometry to hold line segment vertices
          const geometry = new THREE.BufferGeometry();
          // Create a flattened array of vertex coordinates for the line segments
          const vertices = new Float32Array(linesObj.vxs.flatMap(d => [d.x, d.y, d.z]));

          const line = new MeshLine();
          line.setPoints(vertices);

          const material = new MeshLineMaterial({
            color: new THREE.Color(linesObj.color),
            lineWidth: linesObj.displaySize,
            sizeAttenuation: false,  // displaySize is interpreted as pixels.
            resolution: new THREE.Vector2(window.innerWidth, window.innerHeight),
          });

          const mesh = new THREE.Mesh(line, material);
          scene.add(mesh);
        });

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
