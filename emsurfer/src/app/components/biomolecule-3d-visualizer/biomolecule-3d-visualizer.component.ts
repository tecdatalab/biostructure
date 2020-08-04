import { Component, OnInit } from "@angular/core";
import * as NGL from "ngl/dist/ngl.esm.js";

@Component({
  selector: "app-biomolecule-3d-visualizer",
  templateUrl: "./biomolecule-3d-visualizer.component.html",
  styleUrls: ["./biomolecule-3d-visualizer.component.css"],
})
export class BiomoleculeVizualizerComponent implements OnInit {
  stage = null;
  viewer = null;
  constructor() { }

  ngOnInit() {
    this.stage = new NGL.Stage("viewport", {
      theme: "dark",
    });

    this.stage.loadFile(
      "../../../assets/example-files/175d/175d.mrc",
      { defaultRepresentation: true }
    ).then(function (object) {
      object.matrix.elements = [2, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1];

      object.scale.x = 0.5;
      object.scale.y = 0.5;
      object.scale.z = 0.5;

    });

    var shape = new NGL.Shape("shape");
    var sphereBuffer = new NGL.SphereBuffer({
      position: new Float32Array([
        0, 0, 0,
        1, 0, 0,
        2, 0, 0,
        3, 0, 0,
        4, 0, 0,
        5, 0, 0,
        6, 0, 0,
        7, 0, 0,
        8, 0, 0,
        9.5, 9.5, 9.5,

        0, 0, 0,
        0, 1, 0,
        0, 2, 0,
        0, 3, 0,
        0, 4, 0,
        0, 5, 0,
        0, 6, 0,
        0, 7, 0,
        0, 8, 0,
        -25.897052764892578, -39.337318420410156, -25.045808792114258,

        0, 0, 0,
        0, 0, 1,
        0, 0, 2,
        0, 0, 3,
        0, 0, 4,
        0, 0, 5,
        0, 0, 6,
        0, 0, 7,
        0, 0, 8,
        0, 0, 9]),
      color: new Float32Array([
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,
        1, 0, 0,

        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,
        0, 1, 0,

        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1,
        0, 0, 1]),
      radius: new Float32Array([
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,

        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,

        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5])
    });
    shape.addBuffer(sphereBuffer);
    var shapeComp = this.stage.addComponentFromObject(shape);
    shapeComp.addRepresentation("buffer");
    shapeComp.autoView();



    



    this.stage.signals.clicked.add(this.click);

    console.log(this.stage);




  }

  private click(PickingProxy: NGL.pickingProxy) {

    console.log("-----------------");
    console.log(PickingProxy.stage.compList[0].shape.boundingBox);
    console.log("-----------------");




    // console.log(
    //   PickingProxy.picker.surface.center,
    //   PickingProxy.stage.compList[0].position
    // // );
    // PickingProxy.picker.surface.center.x = 0;
    // PickingProxy.picker.surface.center.y = 0;
    // PickingProxy.picker.surface.center.z = 0;

    // if (PickingProxy.mouse.shiftKey) {
    //   const pos = document.getElementById("text") as HTMLElement;
    //   pos.textContent = `x:${PickingProxy.mouse.canvasPosition.x}, y:${PickingProxy.mouse.canvasPosition.y}`;

    //   const shape = new NGL.Shape("shape", { disableImpostor: true });
    //   shape.addSphere(
    //     Array.from({ length: 3 }, () => Math.floor(Math.random() * 600 - 300)),
    //     Array.from({ length: 3 }, () => Math.floor(Math.random() * 250)),
    //     30,
    //     "la bolita"
    //   );
    //   const nosequesesto = PickingProxy.stage.addComponentFromObject(shape);
    //   nosequesesto.addRepresentation("buffer");
    // }

    var centers = function (object: NGL.Viewer) {
      var max = object.boundingBox.max;
      var min = object.boundingBox.min;

      var x_prom = (max.x - min.x) / 2;
      var y_prom = (max.y - min.y) / 2;
      var z_prom = (max.z - min.z) / 2;

      object.boundingBox.max = { x: max.x + x_prom, y: max.y + y_prom, z: max.z + z_prom };
      object.boundingBox.min = { x: min.x + x_prom, y: min.y + y_prom, z: min.z + z_prom };
    }

    centers(PickingProxy.stage.viewer);


    // const shape = new NGL.Shape("shape");
    // var sphereBuffer = new NGL.SphereBuffer({
    //   position: new Float32Array([0, 0, 0]),
    //   color: new Float32Array([255, 0, 0]),
    //   radius: new Float32Array([30])
    // });
    // shape.addBuffer(sphereBuffer);
    // const nosequesesto = PickingProxy.stage.addComponentFromObject(shape);
    // nosequesesto.addRepresentation("buffer");
    // nosequesesto.autoView();

    console.log("-----------------");
    console.log(PickingProxy.stage.compList[0].shape.boundingBox);
    console.log("-----------------");

  }


}
