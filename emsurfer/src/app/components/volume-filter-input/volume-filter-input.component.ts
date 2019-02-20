import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-volume-filter-input',
  templateUrl: './volume-filter-input.component.html',
  styleUrls: ['./volume-filter-input.component.css']
})
export class VolumeFilterInputComponent {

  @Input() parentForm: FormGroup;

  constructor() { }

}
