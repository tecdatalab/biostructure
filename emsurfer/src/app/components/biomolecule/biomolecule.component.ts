import { Component, Input } from "@angular/core";
import { BiomoleculeComparison } from "src/app/models/biomolecule-comparison";
import { Router } from "@angular/router";
import { StringPadder } from 'src/app/models/string-padder';

@Component({
  selector: "app-biomolecule",
  templateUrl: "./biomolecule.component.html",
  styleUrls: ["biomolecule.component.css"]
})
export class BiomoleculeComponent {
  @Input() biomoleculeComparison: BiomoleculeComparison;
  stringPadder: StringPadder;
  constructor(private router: Router) {
    this.stringPadder = new StringPadder();
  }

  searchBiomolecule() {
    const url = "result/" + this.biomoleculeComparison.biomolecule.id;
    this.router.navigate([url], {
      queryParams: { filename: null, mapId: null, contourLevel: null },
      queryParamsHandling: "merge",
      replaceUrl: true
    });
  }
}
