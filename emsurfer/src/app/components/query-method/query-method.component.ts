import { Component, Input, OnInit } from '@angular/core';
import { FormGroup, Validators } from '@angular/forms';

@Component({
  selector: 'app-query-method',
  templateUrl: './query-method.component.html'
})
export class QueryMethodComponent implements OnInit {
  @Input() parentForm: FormGroup;
  cbEmdb: boolean;
  constructor() {}

  ngOnInit() {
    this.cbEmdb = true;
    this.parentForm.get('search_by_emdb_id').valueChanges.subscribe(value => {
      this.cbEmdb = value;
    });
  }

  cbEmdbChange() {
    this.cbEmdb = !this.cbEmdb;
    if (this.cbEmdb) {
      this.parentForm
        .get('emdb_id')
        .setValidators([
          Validators.required,
          Validators.minLength(4),
          Validators.pattern('^[0-9]*$')
        ]);
      this.parentForm
        .get('em_map')
        .get('file')
        .setValidators(null);
      this.parentForm
        .get('em_map')
        .get('contour_level')
        .setValidators(null);
    } else {
      this.parentForm.get('emdb_id').setValidators(null);
      this.parentForm
        .get('em_map')
        .get('file')
        .setValidators(Validators.required);
      this.parentForm
        .get('em_map')
        .get('contour_level')
        .setValidators([
          Validators.required,
          Validators.pattern('^[0-9]+(\.[0-9]+)?$')
        ]);
    }
    this.parentForm.get('em_map').get('file').patchValue(null);
    this.parentForm.get('em_map').get('filename').patchValue(null);
    this.parentForm.get('emdb_id').updateValueAndValidity();
    this.parentForm.get('em_map').get('contour_level').updateValueAndValidity();
    this.parentForm
      .get('em_map')
      .get('file')
      .updateValueAndValidity();
    this.parentForm.patchValue({
      search_by_emdb_id: this.cbEmdb
    });
  }
}
