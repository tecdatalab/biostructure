import { Component, Input } from "@angular/core";
import { BiomoleculeComparison } from "src/app/models/biomolecule-comparison";
import { Router } from "@angular/router";

@Component({
  selector: "app-biomolecule",
  templateUrl: "./biomolecule.component.html",
  styleUrls: ["biomolecule.component.css"]
})
export class BiomoleculeComponent {
  @Input() biomoleculeComparison: BiomoleculeComparison;
  constructor(private router: Router) {}

  searchBiomolecule() {
    const url = "result/" + this.biomoleculeComparison.biomolecule.id;
    this.router.navigate([url], {
      queryParams: { filename: null, mapId: null, contourLevel: null },
      queryParamsHandling: "merge",
      replaceUrl: true
    });
  }
}
