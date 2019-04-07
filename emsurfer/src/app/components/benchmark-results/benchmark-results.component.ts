import { Component, OnInit } from "@angular/core";
import { ActivatedRoute } from "@angular/router";
import { BiomoleculeSearchService } from "src/app/services/biomolecule-search.service";
import { CustomFile } from "src/app/models/custom-file";
import { BenchmarkResult } from "src/app/models/benchmark-result";

@Component({
  selector: "app-benchmark-results",
  templateUrl: "./benchmark-results.component.html"
})
export class BenchmarkResultsComponent implements OnInit {
  files: CustomFile[];
  compressedFilePath: string;
  constructor(
    private route: ActivatedRoute,
    private biomoleculeSearchService: BiomoleculeSearchService
  ) {}

  ngOnInit() {
    const contourRepresentationId = +this.route.snapshot.queryParamMap.get(
      "contourRepresentation"
    );
    const volumeFilter =
      this.route.snapshot.queryParamMap.get("volumeFilter") === "On";
    const topResults = +this.route.snapshot.queryParamMap.get("topResults");
    const emdbIdList = this.route.snapshot.queryParamMap.get("emdbIdList");
    this.biomoleculeSearchService
      .getBatchBiomolecules(
        emdbIdList,
        contourRepresentationId,
        volumeFilter,
        topResults
      )
      .then((data: BenchmarkResult) => {
        this.files = data.results;
        this.files.splice(topResults);
        this.compressedFilePath = data.zipFile;
      });
  }
}
