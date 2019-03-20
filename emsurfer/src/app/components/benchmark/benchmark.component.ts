import { Component, OnInit } from '@angular/core';
import { FormGroup, FormBuilder, Validators } from '@angular/forms';
import { Router } from '@angular/router';
import { FileUploadService } from 'src/app/services/file-upload.service';

@Component({
  selector: 'app-benchmark',
  templateUrl: './benchmark.component.html'
})
export class BenchmarkComponent implements OnInit {
  benchmarkForm: FormGroup;
  defaultFormState: string;
  cbEmdbList: boolean;
  constructor(
    private fb: FormBuilder,
    private router: Router,
    private fileUploadService: FileUploadService
  ) {
    this.cbEmdbList = true;
  }

  ngOnInit() {
    this.benchmarkForm = this.fb.group({
      emdb_id_list: ['1010\n1884\n5502', [Validators.required]],
      file: null,
      contour_representation: 0,
      volume_filter: 'On',
      top_results: 20
    });

    this.defaultFormState = this.benchmarkForm.getRawValue();
  }
  cbEmdbListChange() {
    this.cbEmdbList = !this.cbEmdbList;
    if (this.cbEmdbList) {
      this.benchmarkForm
        .get('emdb_id_list')
        .setValidators([Validators.required]);
      this.benchmarkForm.get('file').setValidators(null);
    } else {
      this.benchmarkForm.get('emdb_id_list').setValidators(null);
      this.benchmarkForm.get('file').setValidators([Validators.required]);
    }
  }

  reset() {
    this.benchmarkForm.reset(this.defaultFormState);
    this.cbEmdbList = true;
  }

  submitHandler() {
    const url = 'benchmark/results';
    const params = {
      emdbIdList: null,
      emdbIdListFile: null,
      contourRepresentation: this.benchmarkForm.get('contour_representation')
        .value,
      volumeFilter: this.benchmarkForm.get('volume_filter').value,
      topResults: this.benchmarkForm.get('top_results').value
    };
    if (this.cbEmdbList) {
      const idListString = this.benchmarkForm.get('emdb_id_list').value;
      params.emdbIdList = this.fileUploadService.uploadIdListString(idListString);
    } else {
      const idListFile = this.benchmarkForm.get('file').value;
      params.emdbIdListFile = this.fileUploadService.uploadIdListFile(idListFile);
    }
    this.router.navigate([url], {
      queryParams: params
    });
  }
}
