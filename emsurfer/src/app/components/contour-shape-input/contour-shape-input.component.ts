import { Component, OnInit, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';
import { ContourRepresentationService } from '../../services/contour-representation.service';

@Component({
  selector: 'app-contour-shape-input',
  templateUrl: './contour-shape-input.component.html'
})
export class ContourShapeInputComponent implements OnInit {
  @Input() parentForm: FormGroup;
  countourRepresentations: Array<string>;
  constructor(
    private contourRepresentationService: ContourRepresentationService
  ) {}

  ngOnInit() {
    this.countourRepresentations = this.contourRepresentationService.getContourShapes();
  }

  dummyFunc() {
    this.countourRepresentations = this.contourRepresentationService.getContourShapes();
  }
}
