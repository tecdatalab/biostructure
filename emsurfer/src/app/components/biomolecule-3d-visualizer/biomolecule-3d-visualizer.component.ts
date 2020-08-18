import { Component, OnInit } from "@angular/core";
import { BiomoleculeVisualizerService } from "../../services/biomolecule-visualizer.service"
import * as NGL from "ngl/dist/ngl.esm.js";
import Color from "olical-color";
import * as THREE from 'three';
import { $ } from 'protractor';

@Component({
  selector: 'app-biomolecule-3d-visualizer',
  templateUrl: './biomolecule-3d-visualizer.component.html',
  styleUrls: ['./biomolecule-3d-visualizer.component.css']
})

export class BiomoleculeVizualizerComponent implements OnInit {

  stage = null;
  structureToShow = null;
  radius = 1;
  segment_file_selected = null;
  points_file_selected = null;
  manualPoints = [
    {
      name: "Group_1", color: "red", points: [
        { name: "1", x: 0.0, y: 0.0, z: 0.0 },
      ]
    },
    {
      name: "Group_2", color: "green", points: [
        { name: "1", x: 0.0, y: 0.0, z: 0.0 }
      ]
    },
    {
      name: "Group_3", color: "blue", points: [
        { name: "1", x: 0.0, y: 0.0, z: 0.0 },
      ]
    }
  ]

  constructor(private BiomoleculeVisualizerService: BiomoleculeVisualizerService) { }

  delay(ms: number) {
    return new Promise(resolve => setTimeout(resolve, ms));
  }

  ngOnInit() {
    var optionToLoad = 7; // number of example option to be requested to the back-end
    // Contruction a set of the stage and the visualizer data
    this.stage = new NGL.Stage("viewport", { theme: "dark" });  // stage where NGL will work
    this.stage.signals.clicked.add(this.click); // definition of click Picking Proxy to interact
    this.loadSystemStructure(optionToLoad); // request and load the data to be shown to the back-end
    this.createInitPoints(); // creation of the inital points
    this.confSideNav("dropdown-btn", null); // set the event listener to the top level obj of the sidebar control panel
  }

  protected createInitPoints() {
    this.manualPoints.forEach(group => {
      group.color = Color.toHex(group.color);
      var color = group.color;
      group.points.forEach(point => {
        this.addPointToStage(group.name + "_" + point.name, point.x, point.y, point.z, Color.toArray(color), this.radius); // test add point
      })
    });
  }

  protected createNewPoint(group_no: number) {
    var gr_points_length = this.manualPoints[group_no].points.length + 1;
    var new_point = { name: "" + gr_points_length, x: 0.0, y: 0.0, z: 0.0 };
    var group = this.manualPoints[group_no];
    this.addPointToStage(group.name + "_" + new_point.name, new_point.x, new_point.y, new_point.z, Color.toArray(group.color), this.radius); // test add point
    this.manualPoints[group_no].points.push(new_point);
  }

  protected deletePoint(group_no: number, point_no: number) {
    this.manualPoints[group_no].points.splice(point_no, 1);
    this.stage.removeComponent(this.stage.getComponentsByName("Group_" + (group_no + 1) + "_" + (point_no + 1)).list[0]);
  }

  protected onPointCoordinateChange(group_no: number, point_no: number, coord: number, value: number) {
    try {
      var stage_point = this.stage.getComponentsByName("Group_" + (group_no + 1) + "_" + (point_no + 1)).list[0];
      if (stage_point) {
        var manual_point = this.manualPoints[group_no].points[point_no];
        switch (coord) {
          case 0:
            manual_point.x = value;
            break;
          case 1:
            manual_point.y = value;
            break;
          case 2:
            manual_point.z = value;
            break
        }
        var position_vector = new THREE.Vector3(manual_point.x, manual_point.y, manual_point.z)
        stage_point.setPosition(position_vector);
      } else {
        console.log("Point doesn't exist");
      }
    } catch (error) {
      console.log(group_no + " " + point_no + " " + coord + " " + value);
      console.log(error);
    }
  }

  protected addPointToStage(shape_name: string, x: number, y: number, z: number, RGBcolor: Array<number>, radius = 1) {
    var shape = new NGL.Shape(shape_name);
    var shapeBuffer = new NGL.SphereBuffer({
      position: new Float32Array([x, y, z]),
      color: new Float32Array([RGBcolor[0], RGBcolor[1], RGBcolor[2]]),
      radius: new Float32Array([radius])
    });
    shape.addBuffer(shapeBuffer);
    var shapeComp = this.stage.addComponentFromObject(shape);
    shapeComp.addRepresentation("buffer");
  }

  protected loadSystemStructure(optionToLoad: number) {
    // request for the structures file to be show
    this.BiomoleculeVisualizerService.getTestStructures(optionToLoad)
      .then((structures: Promise<any>) => {
        this.structureToShow = structures;
        // foreach to load and assaing properties
        Object.keys(structures).forEach(item => {
          var structure = structures[item];
          this.stage.loadFile(structure.path)
            .then(function (object) {
              if (structure.file_type == 0) {
                object.addRepresentation("cartoon", { color: structure.color });
              } else {
                object.addRepresentation("surface", { color: structure.color });
              }

              // console.log(object);

              structure.name = object.name;
              structure.iso = object.reprList[0].repr.isolevel;
              structure.color = Color.toHex(structure.color);

              object.autoView();
            });
        })
      });  // call to get the files
    // console.log(this.stage);
  }

  private click(PickingProxy: NGL.pickingProxy) {

    // var vc_boundingBox = PickingProxy.controls.stage.compList[0].object.boundingBox;
    // var vc_center = PickingProxy.controls.stage.compList[0].object.center;
    // var vc_data = PickingProxy.controls.stage.compList[0].object.data;
    // var mouse_canvasPosition = PickingProxy.mouse.canvasPosition;

    // console.log(PickingProxy);

    // var canvasPosition = PickingProxy.position;  // x , y
    // var z = PickingProxy.controls.viewer.camera.far;

    // var shape = new NGL.Shape("shape");
    // shape.addSphere([canvasPosition.x, canvasPosition.y, canvasPosition.z], [1, 0, 0], 10);
    // var shapeComp = PickingProxy.stage.addComponentFromObject(shape);
    // shapeComp.addRepresentation("buffer");
    // shapeComp.autoView();

  }

  protected onIsoLevelChange(name:string, $event, index:number){
    var component = this.stage.getComponentsByName(name).list[0].reprList[0].repr;
    component.setParameters({ isolevel: $event.target.value })
    this.structureToShow[index].iso = $event.target.value;
  }

  protected onColorChange(name:string, $event, index:number){
    var component = this.stage.getComponentsByName(name).list[0].reprList[0];
    component.setColor($event.target.value);
    this.structureToShow[index].color = $event.target.value;
  }

  protected changeVisibility(name: string, $event) {
    var component = this.stage.getComponentsByName(name).list[0];
    var temp_visibility = component.visible;
    component.setVisibility(!temp_visibility);
    if ($event.target.classList.contains("fa-eye") && !component.visible) {
      $event.target.classList.remove("fa-eye");
      $event.target.classList.add("fa-eye-slash");
    } else if ($event.target.classList.contains("fa-eye-slash") && component.visible) {
      $event.target.classList.remove("fa-eye-slash");
      $event.target.classList.add("fa-eye");
    }
  }

  protected getBackgroundColor(item) {
    return item.color;
  }

  protected confSideNav(class_: string, $event) {
    var dropdown = document.getElementsByClassName(class_);
    for (var i = 0; i < dropdown.length; i++) {
      dropdown[i].addEventListener("click", function () {
        this.classList.toggle("active-item");
        var dropdownContent = this.nextElementSibling;
        if (dropdownContent.style.display === "block") {
          dropdownContent.style.display = "none";
        } else {
          dropdownContent.style.display = "block";
        }
      });
    }
  }

  protected onFileSelected(event) {
    if (event !== undefined) {
      var load_btn, del_btn;
      if (event.target.id === "segment_input") {
        this.segment_file_selected = event.target.files[0];
        load_btn = document.getElementsByName("segments_load_btn")[0];
        del_btn = document.getElementsByName("segments_del_btn")[0];
      } if (event.target.id === "points_input") {
        this.points_file_selected = event.target.files[0];
        load_btn = document.getElementsByName("points_load_btn")[0];
        del_btn = document.getElementsByName("points_del_btn")[0];
      }
      if (load_btn !== undefined || del_btn !== undefined) {
        load_btn.disabled = false;
        del_btn.disabled = false;
      }
    }
  }

  protected onLoadCancel(id: string) {
    var input_file = document.getElementById(id);
    input_file.value = null;
    var load_btn, del_btn;
    if (id === "segment_input") {
      this.segment_file_selected = null;
      load_btn = document.getElementsByName("segments_load_btn")[0];
      del_btn = document.getElementsByName("segments_del_btn")[0];
    } if (id === "points_input") {
      this.points_file_selected = null;
      load_btn = document.getElementsByName("points_load_btn")[0];
      del_btn = document.getElementsByName("points_del_btn")[0];
    }
    if (load_btn !== undefined || del_btn !== undefined) {
      load_btn.disabled = true;
      del_btn.disabled = true;
    }
  }

  protected onLoadConfrim(id: string) {
    var input_file = document.getElementById(id);
    input_file.value = null;
    var load_btn, del_btn;
    if (id === "segment_input") {
      load_btn = document.getElementsByName("segments_load_btn")[0];
      del_btn = document.getElementsByName("segments_del_btn")[0];
      // Load new structure give it
      this.stage.removeAllComponents();
      this.createInitPoints();
      this.stage.loadFile(this.segment_file_selected)
        .then(function (object) {
          object.addRepresentation("surface", {});
          object.autoView();
        });
    } if (id === "points_input") {
      load_btn = document.getElementsByName("points_load_btn")[0];
      del_btn = document.getElementsByName("points_del_btn")[0];

    }
    if (load_btn !== undefined || del_btn !== undefined) {
      load_btn.disabled = true;
      del_btn.disabled = true;
    }
  }

  protected exportSegments(){
    console.log("Segments export");
  }

  protected exportPoints(){
    console.log("Points export");
  }

}