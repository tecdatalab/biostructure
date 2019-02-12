import { Component, OnInit } from '@angular/core';
import { Router } from '@angular/router'
import { FormBuilder, FormGroup, Validators } from '@angular/forms';

@Component({
  selector: 'app-search-form',
  templateUrl: './search-form.component.html',
  styleUrls: ['./search-form.component.css']
})
export class SearchFormComponent implements OnInit {

  searchForm: FormGroup;

  constructor(private fb: FormBuilder, private router: Router) { }

  ngOnInit() {
    
    const countour_representations = this.fb.array([
      "1",
      "2",
      "3"
    ])
   
    const query_group = this.fb.group({
      emdb_id: ['1884',[
        Validators.required,
        Validators.minLength(4)
      ]],
      em_map: this.fb.group({
        file: null,
        contour_level: '3.14'
        })
      })

    const rf_group = this.fb.group({
      min: null,
      max: null
    })

    this.searchForm = this.fb.group({
      contour_shape_representation: countour_representations,
      query: query_group,
      volume_filter: "On",
      resolution_filter: rf_group
    })

    this.searchForm.valueChanges.subscribe(console.log)
  }

  submitHandler(){
    let params = "query/" + 
                this.searchForm.get("query").get("emdb_id").value +
                "/" + this.searchForm.get("volume_filter").value +
                "/" + this.searchForm.get("resolution_filter").get("min").value +
                "/" + this.searchForm.get("resolution_filter").get("max").value
    console.log(params)
    this.router.navigateByUrl(params)
  }
}
