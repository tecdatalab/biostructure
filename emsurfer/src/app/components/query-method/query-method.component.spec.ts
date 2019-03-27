import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { ReactiveFormsModule, FormGroup, FormControl } from "@angular/forms";
import { QueryMethodComponent } from "./query-method.component";
import { EmdbIdInputComponent } from "../emdb-id-input/emdb-id-input.component";
import { UploadEmMapComponent } from "../upload-em-map/upload-em-map.component";

describe("QueryMethodComponent", () => {
  let component: QueryMethodComponent;
  let fixture: ComponentFixture<QueryMethodComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [ReactiveFormsModule],
      declarations: [
        QueryMethodComponent,
        EmdbIdInputComponent,
        UploadEmMapComponent
      ]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(QueryMethodComponent);
    component = fixture.componentInstance;
    component = fixture.componentInstance;
    const dummyParentForm = new FormGroup({
      search_by_emdb_id: new FormControl(),
      emdb_id: new FormControl(),
      em_map: new FormGroup({
        filename: new FormControl(),
        file: new FormControl(),
        contour_level: new FormControl()
      })
    });
    component.parentForm = dummyParentForm;
    fixture.detectChanges();
  });

  it("Radio buttons should call cbEmdbChange() method", async(() => {
    spyOn(component, "cbEmdbChange");
    const rbEmdb = fixture.debugElement.nativeElement.querySelector("#emdb");
    rbEmdb.click();
    fixture.whenStable().then(() => {
      expect(component.cbEmdbChange).toHaveBeenCalled();
    });
  }));
});
