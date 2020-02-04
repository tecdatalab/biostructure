import { Component, Input } from '@angular/core';
import { BiomoleculeComparison } from 'src/app/models/biomolecule-comparison';

@Component({
  selector: 'app-biomolecules-table',
  templateUrl: './biomolecules-table.component.html',
  styleUrls: ["biomolecules-table.component.css"]
})
export class BiomoleculesTableComponent {

  @Input() results: BiomoleculeComparison[];
  constructor() { }

}
