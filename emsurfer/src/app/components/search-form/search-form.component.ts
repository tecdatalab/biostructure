import { Component, OnInit } from '@angular/core';
import { Router } from '@angular/router';
import { FormBuilder, FormGroup, Validators } from '@angular/forms';

@Component({
  selector: 'app-search-form',
  templateUrl: './search-form.component.html'
})
export class SearchFormComponent implements OnInit {
  searchForm: FormGroup;
  defaultFormState: any;
  cbEmdb: boolean;

  constructor(private fb: FormBuilder, private router: Router) {
    this.cbEmdb = true;
  }

  ngOnInit() {
    const queryGroup = this.fb.group({
      emdb_id: [
        '1884',
        [
          Validators.required,
          Validators.minLength(4),
          Validators.pattern('^[0-9]*$')
        ]
      ],
      em_map: this.fb.group({
        file: null,
        contour_level: '3.14'
      })
    });

    const rfGroup = this.fb.group({
      min: null,
      max: null
    });

    this.searchForm = this.fb.group({
      contour_representation: 0,
      query: queryGroup,
      volume_filter: 'On',
      resolution_filter: rfGroup
    });

    this.defaultFormState = this.searchForm.getRawValue();
  }

  submitHandler() {
    const params =
      'query/' +
      this.searchForm.get('query').get('emdb_id').value +
      '/' +
      this.searchForm.get('volume_filter').value +
      '/' +
      this.searchForm.get('resolution_filter').get('min').value +
      '/' +
      this.searchForm.get('resolution_filter').get('max').value;
    this.router.navigateByUrl(params);
  }

  reset() {
    this.searchForm.reset(this.defaultFormState);
    this.cbEmdb = true;
  }
}
