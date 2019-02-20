import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-emdb-id-input',
  templateUrl: './emdb-id-input.component.html'
})
export class EmdbIdInputComponent {
  
  @Input() parentForm: FormGroup;

  constructor() { }

}
