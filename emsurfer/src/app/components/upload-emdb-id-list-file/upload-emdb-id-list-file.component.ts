import { Component, Input, ChangeDetectorRef } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-upload-emdb-id-list-file',
  templateUrl: './upload-emdb-id-list-file.component.html'
})
export class UploadEmdbIdListFileComponent {
  @Input() parentForm: FormGroup;

  constructor(private cd: ChangeDetectorRef) {}

  onFileChange(event: any) {
    const reader = new FileReader();
    if (event.target.files && event.target.files.length) {
      const [file] = event.target.files;
      reader.readAsText(file);
      reader.onload = () => {
        this.parentForm.patchValue({
          file: reader.result
        });
        this.cd.markForCheck();
      };
    } else {
      this.parentForm.get('file').setValue(null);
    }
  }
}
