import { Component, OnInit } from "@angular/core";
import { Router } from "@angular/router";
import { FormBuilder, FormGroup, Validators } from "@angular/forms";
import { FileUploadService } from "src/app/services/file-upload.service";
import { CheckerService } from "src/app/services/checker.service";

@Component({
  selector: "app-search-form",
  templateUrl: "./search-form.component.html",
  styleUrls: ["./search-form.component.css"]
})
export class SearchFormComponent implements OnInit {
  searchForm: FormGroup;
  defaultFormState: string;

  constructor(
    private fb: FormBuilder,
    private router: Router,
    private fileUploadService: FileUploadService,
    private checkerService: CheckerService
  ) {}

  ngOnInit() {
    const queryGroup = this.fb.group({
      search_by_emdb_id: true,
      emdb_id: [
        "1884",
        [
          Validators.required,
          Validators.minLength(4),
          Validators.pattern("^[0-9]*$")
        ]
      ],
      em_map: this.fb.group({
        filename: null,
        file: null,
        contour_level: 3.14
      })
    });

    const rfGroup = this.fb.group({
      min: 0,
      max: 0
    });

    this.searchForm = this.fb.group({
      contour_representation: 0,
      query: queryGroup,
      volume_filter: "Off",
      resolution_filter: rfGroup
    });
    this.defaultFormState = this.searchForm.getRawValue();
  }

  submitHandler() {
    if (this.searchForm.get("query").get("search_by_emdb_id").value) {
      const emdbID = this.searchForm.get("query").get("emdb_id").value;
      this.checkerService.checkBiomolecule(emdbID).then((response: number) => {
        if (response) {
          const url = "result/" + emdbID;
          const params = {
            contourRepresentation: this.searchForm.get("contour_representation")
              .value,
            volumeFilter: this.searchForm.get("volume_filter").value,
            minResolution: this.searchForm.get("resolution_filter").get("min")
              .value,
            maxResolution: this.searchForm.get("resolution_filter").get("max")
              .value
          };
          this.router.navigate([url], {
            queryParams: params
          });
        }
      });
    } else {
      const url = "result/emMap";
      const params = {
        filename: this.searchForm
          .get("query")
          .get("em_map")
          .get("filename").value,
        mapId: null,
        contourLevel: this.searchForm
          .get("query")
          .get("em_map")
          .get("contour_level").value,
        contourRepresentation: this.searchForm.get("contour_representation")
          .value,
        volumeFilter: this.searchForm.get("volume_filter").value,
        minResolution: this.searchForm.get("resolution_filter").get("min")
          .value,
        maxResolution: this.searchForm.get("resolution_filter").get("max").value
      };
      const mapFile = this.searchForm
        .get("query")
        .get("em_map")
        .get("file").value;
      this.fileUploadService
        .uploadEmMap(
          mapFile,
          this.searchForm
            .get("query")
            .get("em_map")
            .get("filename").value
        )
        .then((data: number) => {
          params.mapId = data;
          this.router.navigate([url], {
            queryParams: params
          });
        });
    }
  }

  reset() {
    this.searchForm.reset(this.defaultFormState);
    this.searchForm
      .get("query")
      .get("emdb_id")
      .setValidators([
        Validators.required,
        Validators.minLength(4),
        Validators.pattern("^[0-9]*$")
      ]);
    this.searchForm
      .get("query")
      .get("em_map")
      .get("file")
      .setValidators(null);
    this.searchForm
      .get("query")
      .get("em_map")
      .get("contour_level")
      .setValidators(null);
    this.searchForm
      .get("query")
      .get("em_map")
      .get("file")
      .patchValue(null);
    this.searchForm
      .get("query")
      .get("em_map")
      .get("filename")
      .patchValue(null);
    this.searchForm
      .get("query")
      .get("emdb_id")
      .updateValueAndValidity();
    this.searchForm
      .get("query")
      .get("em_map")
      .get("contour_level")
      .updateValueAndValidity();
    this.searchForm
      .get("query")
      .get("em_map")
      .get("file")
      .updateValueAndValidity();
  }
}
