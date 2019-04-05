import {
  Component,
  OnInit,
  ViewChild,
  ElementRef,
  ɵConsole
} from "@angular/core";
import { ActivatedRoute } from "@angular/router";
import { Chart } from "chart.js";
import { BiomoleculeSearchService } from "../../services/biomolecule-search.service";
import { Biomolecule } from "src/app/models/biomolecule";
import { BiomoleculeComparison } from "src/app/models/biomolecule-comparison";
import { FileDownloadService } from "src/app/services/file-download.service";

@Component({
  selector: "app-search-result",
  templateUrl: "./search-result.component.html",
  providers: [BiomoleculeSearchService]
})
export class SearchResultComponent implements OnInit {
  @ViewChild("canvas") canvasElementRef: ElementRef;
  chart: Chart;
  biomolecule: Biomolecule;
  filename: string;
  filename_result: string;
  results: BiomoleculeComparison[];
  volumeFilter: string;
  isSearchById: boolean;
  downloadResultFile: string;
  descriptors = [];
  values = [];

  constructor(
    private biomoleculeSearchService: BiomoleculeSearchService,
    private fileDownloadService: FileDownloadService,
    private route: ActivatedRoute
  ) {}

  ngOnInit() {
    this.biomolecule = new Biomolecule();
    this.route.params.subscribe(params => {
      this.biomolecule.id = params.emdbId;
      this.load();
    });
  }

  private load() {
    const emdbId = +this.route.snapshot.paramMap.get("emdbId") || null;
    const contourRepresentation = +this.route.snapshot.queryParamMap.get(
      "contourRepresentation"
    );
    const minRes = +this.route.snapshot.queryParamMap.get("minRes");
    const maxRes = +this.route.snapshot.queryParamMap.get("maxRes");
    this.volumeFilter = this.route.snapshot.queryParamMap.get("volumeFilter");
    const mapID = +this.route.snapshot.paramMap.get("mapId");
    if (emdbId) {
      this.biomoleculeSearchService
        .getBiomolecule(emdbId)
        .then((response: Biomolecule) => {
          this.biomolecule = response;
        });
      this.biomoleculeSearchService
        .getZernikeDescriptors(emdbId, contourRepresentation)
        .then(response => {
          this.setValues(response);
        });
      this.biomoleculeSearchService
        .getSimilarBioMolecules(
          emdbId,
          contourRepresentation,
          this.volumeFilter === "On",
          minRes,
          maxRes
        )
        .then(response => {
          this.results = response.results;
          this.filename_result = response.path;
        });
      this.isSearchById = true;
    } else {
      this.filename = this.route.snapshot.queryParamMap.get("filename");
      this.biomoleculeSearchService
        .getZernikeDescriptors(emdbId, contourRepresentation)
        .then(response => {
          this.setValues(response);
        });
      this.biomoleculeSearchService
        .getSimilarBioMoleculesByMap(
          mapID,
          contourRepresentation,
          this.volumeFilter === "On",
          minRes,
          maxRes
        )
        .then(response => {
          this.results = response.results;
          this.filename_result = response.path;
        });
    }
  }

  private initChart(context: ElementRef) {
    this.chart = new Chart(context, {
      type: "line",
      data: {
        labels: this.descriptors,
        datasets: [
          {
            data: this.values,
            borderColor: "black",
            fill: false
          }
        ]
      },
      options: {
        legend: {
          display: false
        },
        scales: {
          xAxes: [
            {
              display: true,
              scaleLabel: {
                display: true,
                labelString: "Zernike Descriptor Number",
                fontSize: 24
              }
            }
          ],
          yAxes: [
            {
              display: true,
              scaleLabel: {
                display: true,
                labelString: "Value",
                fontSize: 24
              }
            }
          ]
        }
      }
    });
  }

  private setValues(value) {
    this.values = value;
    this.descriptors = Array.from(
      new Array(this.values.length),
      (val, index) => index + 1
    );
    const context = this.canvasElementRef.nativeElement;
    this.initChart(context);
  }
}
