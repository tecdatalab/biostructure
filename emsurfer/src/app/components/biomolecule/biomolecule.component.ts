import { Component, Input } from '@angular/core';
import { BiomoleculeComparison } from 'src/app/models/biomolecule-comparison';
import { Router, Params, ActivatedRoute } from '@angular/router';

@Component({
  selector: 'app-biomolecule',
  templateUrl: './biomolecule.component.html',
  styleUrls: ['biomolecule.component.css']
})
export class BiomoleculeComponent {
  @Input() biomoleculeComparison: BiomoleculeComparison;
  constructor(private router: Router, private activatedRoute: ActivatedRoute) {}

  searchBiomolecule() {
    const url = 'result/' + this.biomoleculeComparison.biomolecule.emdb_id;
    this.router.navigate([url], {
      queryParams: { filename: null, mapId: null, contourLevel: null },
      queryParamsHandling: 'merge',
      replaceUrl: true
    });
  }
}
