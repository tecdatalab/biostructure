import { Component, OnInit } from '@angular/core';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls'

import { PDBLoader } from 'three/examples/jsm/loaders/PDBLoader';
import { CSS2DObject } from 'three/examples/jsm/renderers/CSS2DRenderer.js';

@Component({
  selector: 'app-biomolecule-3d-visualizer',
  templateUrl: './biomolecule-3d-visualizer.component.html',
  styleUrls: ['./biomolecule-3d-visualizer.component.css']
})

export class BiomoleculeVizualizerComponent implements OnInit {

  constructor() { }

  ngOnInit() {

    var offset = new THREE.Vector3();

    const aspect_multiplier_x = 1.2;
    const aspect_multiplier_y = 70;

    var scene = new THREE.Scene();
    scene.background = new THREE.Color(0x050505);

    var camera = new THREE.PerspectiveCamera(75, 700 / 480, -40, 1000);
    camera.position.z = 1000;
    scene.add(camera);

    var light = new THREE.DirectionalLight(0xffffff, 0.8);
    light.position.set(1, 1, 1);
    scene.add(light);

    var light = new THREE.DirectionalLight(0xffffff, 0.5);
    light.position.set(- 1, - 1, 1);
    scene.add(light);

    var root = new THREE.Group();
    scene.add(root);

    var renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(window.devicePixelRatio);
    var container = document.querySelector("#three_container")
    renderer.setSize(aspect_multiplier_x * container.clientWidth / 2, aspect_multiplier_y * container.clientHeight / 3);
    container.appendChild(renderer.domElement);

    window.addEventListener('resize', function () {
      renderer.setSize(aspect_multiplier_x * container.clientWidth / 2, aspect_multiplier_y * container.clientHeight / 2);
      camera.aspect = container.clientWidth / container.clientHeight;
      camera.updateProjectionMatrix();
    });

    var controls = new OrbitControls(camera, renderer.domElement);
    controls.update();

    var loader = new PDBLoader();

    while (root.children.length > 0) {

      var object = root.children[0];
      object.parent.remove(object);

    }

    loader.load(

      // map to be loaded
      "../../../assets/example-files/1crn/1crn.pdb",
      // "../../../assets/example-files/emd_1884.map",

      // onLoad callback
      // Here the loaded data is assumed to be an object
      function (pdb) {
        // Add the loaded object to the scene

        console.log(pdb);

        // scene.add(obj);

        var geometryAtoms = pdb.geometryAtoms;
        var geometryBonds = pdb.geometryBonds;
        var json = pdb.json;

        var boxGeometry = new THREE.BoxBufferGeometry(1, 1, 1);
        var sphereGeometry = new THREE.IcosahedronBufferGeometry(1, 2);

        geometryAtoms.computeBoundingBox();
        geometryAtoms.boundingBox.getCenter(offset).negate();

        geometryAtoms.translate(offset.x, offset.y, offset.z);
        geometryBonds.translate(offset.x, offset.y, offset.z);

        var positions = geometryAtoms.getAttribute('position');
        var colors = geometryAtoms.getAttribute('color');

        var position = new THREE.Vector3();
        var color = new THREE.Color();

        for (var i = 0; i < positions.count; i++) {

          position.x = positions.getX(i);
          position.y = positions.getY(i);
          position.z = positions.getZ(i);

          color.r = colors.getX(i);
          color.g = colors.getY(i);
          color.b = colors.getZ(i);

          var material = new THREE.MeshPhongMaterial({ color: color });

          var object = new THREE.Mesh(sphereGeometry, material);
          object.position.copy(position);
          object.position.multiplyScalar(75);
          object.scale.multiplyScalar(25);
          root.add(object);

          var atom = json.atoms[i];

          var text = document.createElement('div');
          text.className = 'label';
          text.style.color = 'rgb(' + atom[3][0] + ',' + atom[3][1] + ',' + atom[3][2] + ')';
          text.textContent = atom[4];

          var label = new CSS2DObject(text);
          label.position.copy(object.position);
          root.add(label);

        }

        positions = geometryBonds.getAttribute('position');

        var start = new THREE.Vector3();
        var end = new THREE.Vector3();

        for (var i = 0; i < positions.count; i += 2) {

          start.x = positions.getX(i);
          start.y = positions.getY(i);
          start.z = positions.getZ(i);

          end.x = positions.getX(i + 1);
          end.y = positions.getY(i + 1);
          end.z = positions.getZ(i + 1);

          start.multiplyScalar(75);
          end.multiplyScalar(75);

          var object = new THREE.Mesh(boxGeometry, new THREE.MeshPhongMaterial(0xffffff));
          object.position.copy(start);
          object.position.lerp(end, 0.5);
          object.scale.set(5, 5, start.distanceTo(end));
          object.lookAt(end);
          root.add(object);

        }

        render();








      },

      // onProgress callback
      function (xhr) {
        console.log((xhr.loaded / xhr.total * 100) + '% loaded');
      },

      // onError callback
      function (err) {
        console.error('An error happened');

        console.log(err);
      }

    );

    // view logic
    var update = function () {
      controls.update();
    };

    // draw scene
    var render = function () {
      renderer.render(scene, camera);
    };

    // run game loop (update, render, repeat)
    var GameLoop = function () {
      requestAnimationFrame(GameLoop);

      update();
      render();
    };

    GameLoop();

  }

}