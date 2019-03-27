import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { ReactiveFormsModule, FormGroup, FormControl } from "@angular/forms";
import { ContourShapeInputComponent } from "./contour-shape-input.component";
import { ContourRepresentationService } from "src/app/services/contour-representation.service";

class MockService {
  getContourShapes() {
    return ["test_array"];
  }
}

describe("ContourShapeInputComponent", () => {
  let component: ContourShapeInputComponent;
  let fixture: ComponentFixture<ContourShapeInputComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [ReactiveFormsModule],
      declarations: [ContourShapeInputComponent],
      providers: [
        {
          provide: ContourRepresentationService,
          useClass: MockService
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
      spyOn(crService, "getContourShapes");
      component.dummyFunc();
      expect(crService.getContourShapes).toHaveBeenCalled();
    }
  ));
});
