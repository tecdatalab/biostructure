import { Component, OnInit } from '@angular/core';
import { FormBuilder, FormGroup, Validators } from '@angular/forms';

@Component({
  selector: 'app-search-form',
  templateUrl: './search-form.component.html',
  styleUrls: ['./search-form.component.css']
})
export class SearchFormComponent implements OnInit {

  searchForm: FormGroup;

  constructor(private fb: FormBuilder) { }

  ngOnInit() {
    
    const countour_representations = this.fb.array([
      "1",
      "2",
      "3"
    ])
   
    const query_group = this.fb.group({
      emdb_id: '1884',
      em_map: this.fb.group({
        file: null,
        contour_level: '3.16'
        })
      })

    const rf_group = this.fb.group({
      min: null,
      max: null
    })

    this.searchForm = this.fb.group({
      contour_shape_representation: countour_representations,
      query: query_group,
      volume_filter: true,
      resolution_filter: rf_group
    })
  }

}
