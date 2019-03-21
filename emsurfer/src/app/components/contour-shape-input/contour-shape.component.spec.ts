import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { FormGroup, FormBuilder } from "@angular/forms";
import { ContourShapeInputComponent } from "./contour-shape-input.component";
import { ContourRepresentationService } from "src/app/services/contour-representation.service";
import { Component } from "@angular/core";

describe("SearchFormComponent", () => {
  let component: ContourShapeInputComponent;
  let fixture: ComponentFixture<ContourShapeInputComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ContourShapeInputComponent, TestComponentWrapper],
      providers: [ContourRepresentationService]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(TestComponentWrapper);
    component = fixture.debugElement.children[0].componentInstance;
    fixture.detectChanges();
  });

  it("ngOnInit should call contourRepresentationService.getContourShapes()", inject(
    [ContourRepresentationService],
    (contourRepresentationService: ContourRepresentationService) => {
      spyOn(contourRepresentationService, "getContourShapes");
      component.ngOnInit();
      fixture.whenStable().then(() => {
        expect(
          contourRepresentationService.getContourShapes
        ).toHaveBeenCalled();
      });
    }
  ));
});

@Component({
  selector: "test-cmp",
  template:
    "<app-contour-shape-input [parentForm]=dummyFormGroup></app-contour-shape-input>"
})
class TestComponentWrapper {
  fb = new FormBuilder();
  dummyFormGroup = this.fb.group({
    contour_representation: 0
  });
}
