import { Component, Input } from '@angular/core';
import { FormGroup, Validators } from '@angular/forms';

@Component({
  selector: 'app-query-method',
  templateUrl: './query-method.component.html'
})
export class QueryMethodComponent {

  @Input() parentForm: FormGroup;
  cbEmdb: boolean;
  constructor() {
    this.cbEmdb = true;
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
    } else {
      this.parentForm.get('emdb_id').setValidators(null);
      this.parentForm
        .get('em_map')
        .get('file')
        .setValidators(Validators.required);
    }
    this.parentForm.get('emdb_id').updateValueAndValidity();
    this.parentForm
      .get('em_map')
      .get('file')
      .updateValueAndValidity();

    this.parentForm.patchValue({
      search_by_emdb_id: this.cbEmdb
    });
  }
}
