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
  defaultFormState: any;
  cb_emdb: boolean;

  constructor(private fb: FormBuilder, private router: Router) { 
    this.cb_emdb = true;
  }

  ngOnInit() {
     
    const query_group = this.fb.group({
      emdb_id: ['1884',[
        Validators.required,
        Validators.minLength(4),
        Validators.pattern("^[0-9]*$")
      ]],
      em_map: this.fb.group({
        file: null,
        contour_level: '3.14'
        })
      })

    const rf_group = this.fb.group({
      min: [null, Validators.required],
      max: [null, Validators.required]
    })

    this.searchForm = this.fb.group({
      contour_representation: 0,
      query: query_group,
      volume_filter: "On",
      resolution_filter: rf_group
    })

    this.defaultFormState = this.searchForm.getRawValue();
    console.log(this.searchForm.get('query').get('emdb_id').value)
    this.searchForm.valueChanges.subscribe(console.log)
  }

  submitHandler(){
    let params = "query/" + 
                this.searchForm.get("query").get("emdb_id").value +
                "/" + this.searchForm.get("volume_filter").value +
                "/" + this.searchForm.get("resolution_filter").get("min").value +
                "/" + this.searchForm.get("resolution_filter").get("max").value
    this.router.navigateByUrl(params)
  }

  cbEmdbChange(){
    this.cb_emdb = !this.cb_emdb
    console.log(this.cb_emdb)
    console.log(this.searchForm.get('query').get('emdb_id').value)
    if (this.cb_emdb){
      this.searchForm.get('query').get('emdb_id').setValidators([
        Validators.required,
        Validators.minLength(4),
        Validators.pattern("^[0-9]*$")
      ])
      this.searchForm.patchValue({
        query: {
          em_map: {
            file : null
          }
        }
      })
    }else{
      this.searchForm.get('query').get('emdb_id').setValidators(null)
      this.searchForm.get('query').get('emdb_id').updateValueAndValidity()
    }
  }

  reset(){
    this.searchForm.patchValue(this.defaultFormState);
    this.cb_emdb = true
  }

}


