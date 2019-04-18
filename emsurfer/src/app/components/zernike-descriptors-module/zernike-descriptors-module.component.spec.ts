import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { ZernikeDescriptorsModuleComponent } from './zernike-descriptors-module.component';

describe('ZernikeDescriptorsModuleComponent', () => {
  let component: ZernikeDescriptorsModuleComponent;
  let fixture: ComponentFixture<ZernikeDescriptorsModuleComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ ZernikeDescriptorsModuleComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ZernikeDescriptorsModuleComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
