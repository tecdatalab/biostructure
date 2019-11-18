import { Component, Input } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-emdb-id-list',
  templateUrl: './emdb-id-list.component.html',
  styleUrls: ['./emdb-id-list.component.css'],
})
export class EmdbIdListComponent {
  @Input() parentForm: FormGroup;

  constructor() {}
}
