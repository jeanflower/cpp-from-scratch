import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js'

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

        const addPointsToScene = (triples) => {
          if (!triples) {
            return;
          }
          // Create a geometry to hold points
          const geometry = new THREE.BufferGeometry();
          // make a flattened map of the coordinates
          const points = new Float32Array(triples.flatMap(d => [d.x, d.y, d.z]));
          // console.log(`point data is ${points}`);
          geometry.setAttribute("position", new THREE.BufferAttribute(points, 3));
          // Create a material for the points
          const material = new THREE.PointsMaterial({ color: 0x00ff00, size: 0.03 });
          const pointCloud = new THREE.Points(geometry, material);
          // Add the points to the scene
          scene.add(pointCloud);
        }
        addPointsToScene(data.pts);

        const addPolylineToScene = (triples) => {
          // Create a geometry to hold line segment vertices
          const geometry = new THREE.BufferGeometry();
          // Create a flattened array of vertex coordinates for the line segments
          const vertices = new Float32Array(triples.flatMap(d => [d.x, d.y, d.z]));
          // console.log(`polyline data is ${vertices}`);
          geometry.setAttribute("position", new THREE.BufferAttribute(vertices, 3));
      
          // Create a material for the line segments
          const material = new THREE.LineBasicMaterial({ color: 0x0000ff, linewidth: 2 });
      
          // Create the line segments and add them to the scene
          const lineSegments = new THREE.Line(geometry, material);
          scene.add(lineSegments);
        };
        data.polylines.map(addPolylineToScene);

        // Animate the scene
        animate();
    });

// Set the camera position
camera.position.z = 20;

// Animation loop
function animate() {
    requestAnimationFrame(animate);
    // controls.update();
    renderer.render(scene, camera);
}
