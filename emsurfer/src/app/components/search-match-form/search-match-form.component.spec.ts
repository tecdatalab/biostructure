import { async, ComponentFixture, TestBed } from '@angular/core/testing';

import { SearchMatchFormComponent } from './search-match-form.component';

describe('SearchMatchFormComponent', () => {
  let component: SearchMatchFormComponent;
  let fixture: ComponentFixture<SearchMatchFormComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [ SearchMatchFormComponent ]
    })
    .compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(SearchMatchFormComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
