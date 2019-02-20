import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-resolution-filter',
  templateUrl: './resolution-filter.component.html'
})
export class ResolutionFilterComponent  {

  @Input() parentForm: FormGroup;
  constructor() { }

}
