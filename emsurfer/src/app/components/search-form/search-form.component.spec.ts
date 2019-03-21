import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from '@angular/core/testing';
import { ReactiveFormsModule } from '@angular/forms';
import { RouterTestingModule } from '@angular/router/testing';
import { ContourShapeInputComponent } from '../contour-shape-input/contour-shape-input.component';
import { VolumeFilterInputComponent } from '../volume-filter-input/volume-filter-input.component';
import { SearchFormComponent } from './search-form.component';
import { HeaderComponent } from '../header/header.component';
import { EmdbIdInputComponent } from '../emdb-id-input/emdb-id-input.component';
import { UploadEmMapComponent } from '../upload-em-map/upload-em-map.component';
import { QueryMethodComponent } from '../query-method/query-method.component';
import { ResolutionFilterComponent } from '../resolution-filter/resolution-filter.component';
import { Router } from '@angular/router';

class MockRouter {
  navigate(urls: string[], extras: string) {
    return true;
  }
}

describe('SearchFormComponent', () => {
  let component: SearchFormComponent;
  let fixture: ComponentFixture<SearchFormComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [
        ReactiveFormsModule,
        RouterTestingModule
      ],
      declarations: [
        SearchFormComponent,
        HeaderComponent,
        ContourShapeInputComponent,
        VolumeFilterInputComponent,
        EmdbIdInputComponent,
        UploadEmMapComponent,
        QueryMethodComponent,
        ResolutionFilterComponent,
      ],
      providers: [
        { provide: Router, useClass: MockRouter }
      ]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SearchFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('submit button should call submitHandler()', async(() => {
    spyOn(component, 'submitHandler');
    const button = fixture.debugElement.nativeElement.querySelector('#submitButton');
    button.click();
    fixture.whenStable().then(() => {
      expect(component.submitHandler).toHaveBeenCalled();
    });
  }));

  it('reset button should call reset', async(() => {
    spyOn(component, 'reset');
    const button = fixture.debugElement.nativeElement.querySelector('#resetButton');
    button.click();
    fixture.whenStable().then(() => {
      expect(component.reset).toHaveBeenCalled();
    });
  }));

  it('reset function should call this.searchForm.reset(expectedParameter)', async(() => {
    spyOn(component.searchForm, 'reset');
    component.reset();
    const expectedParameter = component.defaultFormState;
    const expectedUrl = 'result/1884';
    expect(component.searchForm.reset).toHaveBeenCalledWith(expectedParameter);
  }));


  it('submitHandler should call router navigate with expected query params', inject([Router], ( router: Router ) => {
    spyOn(router, 'navigate');
    component.submitHandler();
    const expectedQueryParams = {
      contourRepresentation: 0,
      volumeFilter: 'On',
      minResolution: null,
      maxResolution: null
    };
    const expectedUrl = 'result/1884';
    expect(router.navigate).toHaveBeenCalledWith([expectedUrl], {queryParams: expectedQueryParams} );
  }));

  it('submitHandler should call router navigate with expected query params if the user wants to do a query with a file', 
  inject([Router], ( router: Router ) => {
    spyOn(router, 'navigate');
    component.searchForm.patchValue({
      query:{
        search_by_emdb_id: false,
        em_map:{
          filename: 'test_file.map',
        }
      }
    });
    component.submitHandler();
    const expectedQueryParams = {
      filename: 'test_file.map',
      mapId: 4444,
      contourLevel: 3.14,
      contourRepresentation: 0,
      volumeFilter: 'On',
      minResolution: null,
      maxResolution: null
    };
    const expectedUrl = 'result/emMap';
    expect(router.navigate).toHaveBeenCalledWith([expectedUrl], {queryParams: expectedQueryParams} );
  }));
});
