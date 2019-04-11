import { Component, OnInit, Input } from "@angular/core";
import { FormGroup } from "@angular/forms";
import { ContourRepresentationService } from "../../services/contour-representation.service";
import { ContourRepresentation } from "../../models/contour-representation";

@Component({
  selector: "app-contour-shape-input",
  templateUrl: "./contour-shape-input.component.html"
})
export class ContourShapeInputComponent implements OnInit {
  @Input() parentForm: FormGroup;
  contourRepresentations: ContourRepresentation[];
  constructor(
    private contourRepresentationService: ContourRepresentationService
  ) {}

  ngOnInit() {
    this.contourRepresentationService
      .getContourShapes()
      .then((data: ContourRepresentation[]) => {
        this.contourRepresentations = data;
      });
  }
}
