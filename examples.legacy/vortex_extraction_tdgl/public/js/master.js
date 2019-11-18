var clock = new THREE.Clock();
var scene = new THREE.Scene();

var camera = new THREE.PerspectiveCamera(30, window.innerWidth/window.innerHeight, 0.1, 1000); 
camera.position.z = 200;

var renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
renderer.setPixelRatio(window.devicePixelRatio);
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.setClearColor(0xffffff, 1);
document.body.appendChild(renderer.domElement);

var pointLight = new THREE.PointLight(0xffffff);
pointLight.position.x = 100;
pointLight.position.y = 100;
pointLight.position.z = 100;
// pointLight.castShadow = true;
// pointLight.shadowDarkness = 0.5;
scene.add(pointLight);

var directionalLight = new THREE.DirectionalLight(0xffffff);
scene.add(directionalLight);

cameraControls = new THREE.TrackballControls(camera, renderer.domElement);
cameraControls.target.set(0, 0, 0);
cameraControls.zoomSpeed = 0.04;
cameraControls.panSpeed = 0.04;
// cameraControls.addEventListener("change", render); // not working.. sigh

var raycaster = new THREE.Raycaster();
var mousePos = new THREE.Vector2();

window.addEventListener("mousedown", onMouseDown, false );
window.addEventListener("mousemove", onMouseMove, false);
window.addEventListener("resize", onResize, false);

function onMouseDown(evt) {
  mousePos.x = (evt.clientX / window.innerWidth) * 2 - 1;
  mousePos.y = -(evt.clientY / window.innerHeight) * 2 + 1;
}

function onMouseMove(evt) {
  mousePos.x = (evt.clientX / window.innerWidth) * 2 - 1;
  mousePos.y = -(evt.clientY / window.innerHeight) * 2 + 1;
}

function onResize() {
  camera.aspect = window.innerWidth/window.innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(window.innerWidth, window.innerHeight);
  cameraControls.handleResize();
}

function render() {
  raycaster.setFromCamera(mousePos, camera);
  // var intersects = raycaster.intersectObjects(scene.children);

  // scene
  var delta = clock.getDelta();
  requestAnimationFrame(render);
  cameraControls.update(delta);
  directionalLight.position.copy(camera.position);
  renderer.render(scene, camera);
}

const radius = 0.2;

$.getJSON("vortices", function(data) {
  vortices = data.vortices;
  for (i = 0; i < vortices.length; i ++) {
    var verts = vortices[i];
    console.log(verts);

    // var lineGeometry = new THREE.Geometry();
    var points = [];
    for (j = 0; j < verts.length; j ++) {
      points.push(new THREE.Vector3(verts[j].x, verts[j].y, verts[j].z));
      // lineGeometry.vertices.push(new THREE.Vector3(verts[j].x, verts[j].y, verts[j].z));
    }
    var curve = new THREE.CatmullRomCurve3(points);

    // var lineMaterial = new THREE.LineBasicMaterial({color: 'black'});
    // var line = new THREE.Line(lineGeometry, lineMaterial);
    // scene.add(line);

    var tubeGeometry = new THREE.TubeGeometry(curve, 100, radius, 8, false);
    var tubeMaterial = new THREE.MeshPhysicalMaterial({
      color: 'black',
      side: THREE.DoubleSide,
      wireframe: false
    });
    var tubeMesh = new THREE.Mesh(tubeGeometry, tubeMaterial);
    scene.add(tubeMesh);
  }

  render();
});
  
