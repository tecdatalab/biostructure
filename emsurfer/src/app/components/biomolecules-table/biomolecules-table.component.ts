import { Component, Input } from '@angular/core';
import { BiomoleculeComparison } from 'src/app/models/biomolecule-comparison';

@Component({
  selector: 'app-biomolecules-table',
  templateUrl: './biomolecules-table.component.html'
})
export class BiomoleculesTableComponent {

  @Input() results: BiomoleculeComparison[];
  constructor() { }

}
