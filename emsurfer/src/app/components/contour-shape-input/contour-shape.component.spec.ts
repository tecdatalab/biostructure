import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { ReactiveFormsModule, FormGroup, FormControl } from "@angular/forms";
import { ContourShapeInputComponent } from "./contour-shape-input.component";
import { ContourRepresentationService } from "src/app/services/contour-representation.service";
import { ContourRepresentation } from "src/app/models/contour-representation";
import { HttpClientModule } from "@angular/common/http";

class MockPromise {
  then(func) {
    const contourRepresentation = new ContourRepresentation();
    contourRepresentation.id = 0;
    contourRepresentation.name = "test";
    func([contourRepresentation]);
  }
}

class MockContourRepresentationService {
  getContourShapes() {
    return new MockPromise();
  }
}

describe("ContourShapeInputComponent", () => {
  let component: ContourShapeInputComponent;
  let fixture: ComponentFixture<ContourShapeInputComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [ReactiveFormsModule, HttpClientModule],
      declarations: [ContourShapeInputComponent],
      providers: [
        {
          provide: ContourRepresentationService,
          useClass: MockContourRepresentationService
        }
      ]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ContourShapeInputComponent);
    component = fixture.componentInstance;
    const dummyParentForm = new FormGroup({
      contour_representation: new FormControl()
    });
    component.parentForm = dummyParentForm;
    fixture.detectChanges();
  });

  it("ngOnInit should call contourRepresentationService.getContourShapes()", inject(
    [ContourRepresentationService],
    (crService: ContourRepresentationService) => {
      spyOn(crService, "getContourShapes").and.returnValue(new MockPromise());
      component.ngOnInit();
      expect(crService.getContourShapes).toHaveBeenCalled();
    }
  ));
});
