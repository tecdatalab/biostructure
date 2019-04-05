import { Component, Input, ChangeDetectorRef } from '@angular/core';
import { FormGroup } from '@angular/forms';

@Component({
  selector: 'app-upload-em-map',
  templateUrl: './upload-em-map.component.html'
})
export class UploadEmMapComponent {
  @Input() parentForm: FormGroup;

  constructor(private cd: ChangeDetectorRef) {}

  onFileChange(event: any) {
    const reader = new FileReader();
    if (event.target.files && event.target.files.length) {
      const [file] = event.target.files;
      reader.readAsDataURL(file);
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
