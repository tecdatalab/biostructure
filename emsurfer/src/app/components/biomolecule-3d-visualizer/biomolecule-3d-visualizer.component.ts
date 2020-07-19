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
  constructor() {}

  ngOnInit() {
    this.stage = new NGL.Stage("viewport", {
      theme: "dark",
    });

    let schemeId = NGL.ColormakerRegistry.addSelectionScheme(
      [
        ["red", "64-74 or 134-154 or 222-254 or 310-310 or 322-326"],
        ["green", "311-322"],
        [
          "yellow",
          "40-63 or 75-95 or 112-133 or 155-173 or 202-221 or 255-277 or 289-309",
        ],
        ["blue", "1-39 or 96-112 or 174-201 or 278-288"],
        ["white", "*"],
      ],
      "Transmembrane 3dqb"
    );

    this.stage.loadFile(
      //   "http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/map/emd_0001.map.gz",
      // "../../../assets/example-files/emd_1884.map",
      //   "../../../assets/example-files/emd_0001.map.gz",
      "../../../assets/example-files/emd_1884.map",
      { defaultRepresentation: true }
    );

    this.stage.signals.clicked.add(this.click);
    console.log(this.stage);
  }

  private click(PickingProxy: NGL.pickingProxy) {
    console.log(
      PickingProxy.picker.surface.center,
      PickingProxy.stage.compList[0].position
    );
    //PickingProxy.picker.surface.center.x = 0;
    //PickingProxy.picker.surface.center.y = 0;
    //PickingProxy.picker.surface.center.z = 0;

    if (PickingProxy.mouse.shiftKey) {
      const pos = document.getElementById("text") as HTMLElement;
      pos.textContent = `x:${PickingProxy.mouse.canvasPosition.x}, y:${PickingProxy.mouse.canvasPosition.y}`;

      const shape = new NGL.Shape("shape", { disableImpostor: true });
      shape.addSphere(
        Array.from({ length: 3 }, () => Math.floor(Math.random() * 600 - 300)),
        Array.from({ length: 3 }, () => Math.floor(Math.random() * 250)),
        30,
        "la bolita"
      );
      const nosequesesto = PickingProxy.stage.addComponentFromObject(shape);
      nosequesesto.addRepresentation("buffer");
    }
  }
}
