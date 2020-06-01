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
    this.stage.loadFile(
      "http://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-0001/map/emd_0001.map.gz",
      { defaultRepresentation: true }
    );
  }
}
