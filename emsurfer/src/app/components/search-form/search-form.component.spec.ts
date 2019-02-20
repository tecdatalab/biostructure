import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from '@angular/core/testing';
import { ReactiveFormsModule } from '@angular/forms';
import { RouterTestingModule } from '@angular/router/testing';
import { ContourShapeInputComponent } from '../contour-shape-input/contour-shape-input.component';
import { SearchFormComponent } from './search-form.component';
import { HeaderComponent } from '../header/header.component';
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
  let de: DebugElement;

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
        DummyComponent
      ]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SearchFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('simulate button click: should redirect to default url on submit', async(
    inject([Router, Location], (router: Router, location: Location) => {
      fixture.debugElement.query(By.css('button')).nativeElement.click();
      fixture.whenStable().then(() => {
        expect(location.path()).toEqual('/query/1884/On/' + null + '/' + null);
      });
    })
  ));
});
