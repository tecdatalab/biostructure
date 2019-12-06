import { Component, OnInit, Input } from "@angular/core";
import { FormGroup } from "@angular/forms";
import { ContourRepresentationService } from "../../services/contour-representation.service";
import { ContourRepresentation } from "../../models/contour-representation";

@Component({
  selector: "app-contour-shape-input",
  templateUrl: "./contour-shape-input.component.html",
  styleUrls: ['./contour-shape-input.component.css'],
})
export class ContourShapeInputComponent implements OnInit {
  @Input() parentForm: FormGroup;
  contourRepresentations = [{ id: 1 }];
  constructor(
    private contourRepresentationService: ContourRepresentationService
  ) {}

  ngOnInit() {
    sessionStorage.typeDescriptor = 1;
    this.contourRepresentationService
      .getContourShapes()
      .then((data: ContourRepresentation[]) => {
        this.contourRepresentations = data;
      });
  }

  onChange(value){
    sessionStorage.typeDescriptor = value[0];
  }
  
}
