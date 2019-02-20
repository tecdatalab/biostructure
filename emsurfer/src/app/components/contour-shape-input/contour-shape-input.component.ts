import { Component, OnInit, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-contour-shape-input',
  templateUrl: './contour-shape-input.component.html',
  styleUrls: ['./contour-shape-input.component.css']
})
export class ContourShapeInputComponent implements OnInit {
  @Input() parentForm: FormGroup;
  countourRepresentations: Array<string>;
  constructor() {}

  ngOnInit() {
    this.countourRepresentations = [
      'EMDB contour',
      'EMDB contour + 1/3 core',
      'EMDB contour + 2/3 core',
      'EMDB contour + 1/3 + 2/3 core',
      'EMDB contour + 1 std dev'
    ];
  }
}
