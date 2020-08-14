import { Component, OnInit } from "@angular/core";
import { BiomoleculeVisualizerService } from "../../services/biomolecule-visualizer.service"
import * as NGL from "ngl/dist/ngl.esm.js";
import { importExpr } from '@angular/compiler/src/output/output_ast';
import { reduce } from 'rxjs/operators';

@Component({
  selector: 'app-biomolecule-3d-visualizer',
  templateUrl: './biomolecule-3d-visualizer.component.html',
  styleUrls: ['./biomolecule-3d-visualizer.component.css']
})

export class BiomoleculeVizualizerComponent implements OnInit {

  stage = null;
  structureToShow = null;
  manualPoints = [
    {
      name: "Group_1", color: "red", points: [
        { x: 0.0, y: 0.0, z: 0.0 },
        { x: 1.0, y: 1.0, z: 1.0 },
      ]
    },
    {
      name: "Group_2", color: "red", points: [
        { x: 0.0, y: 0.0, z: 0.0 }
      ]
    },
  ]

  constructor(private BiomoleculeVisualizerService: BiomoleculeVisualizerService) { }

  delay(ms: number) {
    return new Promise(resolve => setTimeout(resolve, ms));
  }

  ngOnInit() {

    var temp_delay = 400;   // delay time between loading files
    this.stage = new NGL.Stage("viewport", { theme: "dark" });    // stage where NGL will work
    this.stage.signals.clicked.add(this.click);   // definition of click Picking Proxy to interact
    // request for the structures file to be show
    var optionToLoad = 7;
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
              object.autoView();
            });
        })
      });  // call to get the files
    this.confSideNav("dropdown-btn", null);
  }

  private click(PickingProxy: NGL.pickingProxy) {
    var vc_boundingBox = PickingProxy.controls.stage.compList[0].object.boundingBox;
    var vc_center = PickingProxy.controls.stage.compList[0].object.center;
    var vc_data = PickingProxy.controls.stage.compList[0].object.data;
    var mouse_canvasPosition = PickingProxy.mouse.canvasPosition;

    console.log(PickingProxy);

    var canvasPosition = PickingProxy.position;  // x , y
    var z = PickingProxy.controls.viewer.camera.far;

    var shape = new NGL.Shape("shape");
    shape.addSphere([canvasPosition.x, canvasPosition.y, canvasPosition.z], [1, 0, 0], 10);
    var shapeComp = PickingProxy.stage.addComponentFromObject(shape);
    shapeComp.addRepresentation("buffer");
    shapeComp.autoView();

  }

  getBackgroundColor(item) {
    return item.color;
  }

  confSideNav(class_: string, $event) {
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


}