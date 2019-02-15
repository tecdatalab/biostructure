import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { VolumeFilterInputComponent } from './volume-filter-input.component';

describe('VolumeFilterInputComponent', () => {
  let component: VolumeFilterInputComponent;
  let fixture: ComponentFixture<VolumeFilterInputComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ VolumeFilterInputComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(VolumeFilterInputComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

});
