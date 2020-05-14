import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-resolution-filter',
  templateUrl: './resolution-filter.component.html',
  styleUrls: ['./resolution-filter.component.css'],
})
export class ResolutionFilterComponent  {

  @Input() parentForm: FormGroup;
  constructor() { }

}
