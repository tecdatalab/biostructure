import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { ZernikeResultComponent } from './zernike-result.component';

describe('ZernikeResultComponent', () => {
  let component: ZernikeResultComponent;
  let fixture: ComponentFixture<ZernikeResultComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ ZernikeResultComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ZernikeResultComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
