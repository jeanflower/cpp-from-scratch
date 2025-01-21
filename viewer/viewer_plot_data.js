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
fetch("output/xyz_coordinates.json")
    .then(response => 
      {
        console.log(`got response from xyz_coordinates.json file ${JSON.stringify(response)}`)
        return response.json();
      }
    )
    .then(data => {
        // Create a geometry to hold points
        const geometry = new THREE.BufferGeometry();
        const points = new Float32Array(data.flatMap(d => [d.x, d.y, d.z]));
        geometry.setAttribute("position", new THREE.BufferAttribute(points, 3));

        // Create a material for the points
        const material = new THREE.PointsMaterial({ color: 0x00ff00, size: 0.05 });
        const pointCloud = new THREE.Points(geometry, material);

        // Add the points to the scene
        scene.add(pointCloud);
                
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
