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
import { DebugElement } from '@angular/core';
import { By } from '@angular/platform-browser';
import { Location, CommonModule } from '@angular/common';
import { Component } from '@angular/core';
import { Router } from '@angular/router';

@Component({
  template: ''
})
class DummyComponent {}

describe('SearchFormComponent', () => {
  let component: SearchFormComponent;
  let fixture: ComponentFixture<SearchFormComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      imports: [
        ReactiveFormsModule,
        RouterTestingModule,
        CommonModule,
        RouterTestingModule.withRoutes([
          { path: 'query/:id/:volume/:min/:max', component: DummyComponent }
        ])
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
        DummyComponent
      ]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SearchFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('check validators of searchform', async(() => {
    component.searchForm.get('query').get('emdb_id').setValue(null);
    expect(component.searchForm.valid).toBeFalsy();
    component.searchForm.get('query').get('emdb_id').setValue('188');
    expect(component.searchForm.valid).toBeFalsy();
    component.searchForm.get('query').get('emdb_id').setValue('188A');
    expect(component.searchForm.valid).toBeFalsy();
  }));

  it('simulate button click: should redirect to query url on submit', async(
    inject([Router, Location], (router: Router, location: Location) => {
      fixture.debugElement.query(By.css('button')).nativeElement.click();
      fixture.whenStable().then(() => {
        const query = component.searchForm.get('query').get('emdb_id').value;
        const filter = component.searchForm.get('volume_filter').value;
        expect(location.path()).toEqual('/query/' + query + '/' + filter + '/' + null + '/' + null);
      });
    })
  ));
})