import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-volume-filter-input',
  templateUrl: './volume-filter-input.component.html'
})
export class VolumeFilterInputComponent {

  @Input() parentForm: FormGroup;

  constructor() { }

}
