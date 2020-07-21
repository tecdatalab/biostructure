import { Component, OnInit } from "@angular/core";
import { BiomoleculeVisualizerService } from "../../services/biomolecule-visualizer.service"
import * as NGL from "ngl/dist/ngl.esm.js";
import { importExpr } from '@angular/compiler/src/output/output_ast';

@Component({
  selector: 'app-biomolecule-3d-visualizer',
  templateUrl: './biomolecule-3d-visualizer.component.html',
  styleUrls: ['./biomolecule-3d-visualizer.component.css']
})

export class BiomoleculeVizualizerComponent implements OnInit {

  stage = null;
  structuresToLoads = null;
  structureToShow = [];

  constructor(private BiomoleculeVisualizerService: BiomoleculeVisualizerService) { 
    // request for the structures file to be show
    var optionToLoad = 7;
    this.structuresToLoads = this.BiomoleculeVisualizerService.getTestStructures(optionToLoad);  // call to get the files
    Object.keys(this.structuresToLoads).forEach(item => {
      this.structureToShow.push(this.structuresToLoads[item]);
    });
    console.log(this.structureToShow);
  }

  delay(ms: number) {
    return new Promise(resolve => setTimeout(resolve, ms));
  }

  ngOnInit() {

    var temp_delay = 400;   // delay time between loading files
    this.stage = new NGL.Stage("viewport", { theme: "dark" });    // stage where NGL will work
    this.stage.signals.clicked.add(this.click);   // definition of click Picking Proxy to interact
    
    // foreach to load and assaing properties
    Object.keys(this.structuresToLoads).forEach(item => {
      var structure = this.structuresToLoads[item];
      this.stage.loadFile(structure.path)
        .then(function (object) {
          if (structure.file_type == 0) {
            object.addRepresentation("cartoon", { color: structure.color });
          } else {
            object.addRepresentation("surface", { color: structure.color });
          }
          object.autoView();
        });
      this.delay(temp_delay)
    })



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
    shape.addSphere([canvasPosition.x, canvasPosition.y, canvasPosition.z], [ 1, 0, 0 ], 10);
    var shapeComp = PickingProxy.stage.addComponentFromObject(shape);
    shapeComp.addRepresentation("buffer");
    shapeComp.autoView();

  }


  getBackgroundColor(item){
    return item.color;
  }
}
