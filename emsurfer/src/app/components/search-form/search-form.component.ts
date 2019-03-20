import { Component, OnInit } from '@angular/core';
import { Router, NavigationExtras } from '@angular/router';
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
      search_by_emdb_id: true,
      emdb_id: [
        '1884',
        [
          Validators.required,
          Validators.minLength(4),
          Validators.pattern('^[0-9]*$')
        ]
      ],
      em_map: this.fb.group({
        filename: null,
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
    if (this.searchForm.get('query').get('search_by_emdb_id').value) {
      const url = 'result/' +  this.searchForm.get('query').get('emdb_id').value;
      const params = {
        contourRepresentation: this.searchForm.get('contour_representation')
          .value,
        volumeFilter: this.searchForm.get('volume_filter').value,
        minResolution: this.searchForm.get('resolution_filter').get('min')
          .value,
        maxResolution: this.searchForm.get('resolution_filter').get('max').value
      };
      this.router.navigate([url], {
        queryParams: params
      });
    } else {
      const url =
        'result/emMap';
      const params = {
        filename: this.searchForm
          .get('query')
          .get('em_map')
          .get('filename').value,
        contourLevel: this.searchForm
          .get('query')
          .get('em_map')
          .get('contour_level').value,
        contourRepresentation: this.searchForm.get('contour_representation')
          .value,
        volumeFilter: this.searchForm.get('volume_filter').value,
        minResolution: this.searchForm.get('resolution_filter').get('min')
          .value,
        maxResolution: this.searchForm.get('resolution_filter').get('max').value
      };
      this.router.navigate([url], {
        queryParams: params
      });
    }
  }

  reset() {
    this.searchForm.reset(this.defaultFormState);
    this.cbEmdb = true;
  }
}
