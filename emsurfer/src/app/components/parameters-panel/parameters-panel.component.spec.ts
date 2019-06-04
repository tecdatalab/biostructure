import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { ParametersPanelComponent } from './parameters-panel.component';

describe('ParametersPanelComponent', () => {
  let component: ParametersPanelComponent;
  let fixture: ComponentFixture<ParametersPanelComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ ParametersPanelComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(ParametersPanelComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });
  /*
  it('should create', () => {
    expect(component).toBeTruthy();
  });*/
});
