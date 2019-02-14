import { Component, OnInit, Input } from '@angular/core';
import { FormGroup, FormControl } from '@angular/forms';

@Component({
  selector: 'app-volume-filter-input',
  templateUrl: './volume-filter-input.component.html',
  styleUrls: ['./volume-filter-input.component.css']
})
export class VolumeFilterInputComponent implements OnInit {

  @Input() parentForm: FormGroup;

  constructor() { }

  ngOnInit() {
  }

}
