import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { ContourShapeInputComponent } from './contour-shape-input.component';

describe('ContourShapeInputComponent', () => {
  let component: ContourShapeInputComponent;
  let fixture: ComponentFixture<ContourShapeInputComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ ContourShapeInputComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ContourShapeInputComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
