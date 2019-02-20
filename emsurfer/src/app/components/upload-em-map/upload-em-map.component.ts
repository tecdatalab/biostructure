import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-upload-em-map',
  templateUrl: './upload-em-map.component.html'
})
export class UploadEmMapComponent {

  @Input() parentForm: FormGroup;

  constructor() { }


}
