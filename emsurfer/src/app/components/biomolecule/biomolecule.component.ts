import { Component, Input } from '@angular/core';
import { BiomoleculeComparison } from 'src/app/models/biomolecule-comparison';

@Component({
  selector: 'app-biomolecule',
  templateUrl: './biomolecule.component.html'
})
export class BiomoleculeComponent  {

  @Input() biomoleculeComparison: BiomoleculeComparison;
  constructor() { }

}
