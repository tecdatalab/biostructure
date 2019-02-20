import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

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

  cbEmdbChange(){
    this.cbEmdb = !this.cbEmdb;
  }
}
