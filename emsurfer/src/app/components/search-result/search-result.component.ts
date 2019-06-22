import { Component, OnInit, ViewChild, ElementRef } from "@angular/core";
import { ActivatedRoute } from "@angular/router";
import { Chart } from "chart.js";
import { BiomoleculeSearchService } from "../../services/biomolecule-search.service";
import { Biomolecule } from "src/app/models/biomolecule";
import { BiomoleculeComparison } from "src/app/models/biomolecule-comparison";
import { DescriptorService } from "src/app/services/descriptor.service";
import { StringPadder } from "src/app/models/string-padder";
import { Descriptor } from "src/app/models/descriptor";

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
  volumeFilter: string;
  results: BiomoleculeComparison[];
  isSearchById: boolean;
  downloadResultFile: string;
  descriptors = [];
  values = [];
  stringPadder: StringPadder;

  constructor(
    private biomoleculeSearchService: BiomoleculeSearchService,
    private descriptorService: DescriptorService,
    private route: ActivatedRoute
  ) {
    this.stringPadder = new StringPadder();
  }

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
    const minRes = this.route.snapshot.queryParamMap.get("minResolution");
    const maxRes = this.route.snapshot.queryParamMap.get("maxResolution");
    this.volumeFilter = this.route.snapshot.queryParamMap.get("volumeFilter");
    if (emdbId) {
      this.biomoleculeSearchService
        .getBiomolecule(emdbId)
        .then((biomoleculeResponse: Biomolecule) => {
          if (biomoleculeResponse) {
            this.biomolecule = biomoleculeResponse;
            this.descriptorService
              .getDescriptor(emdbId, contourRepresentation)
              .then((response: Descriptor) => {
                this.setValues(response.numbers);
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
                this.downloadResultFile = response.path;
              });
          }
        });
      this.isSearchById = true;
    } else {
      //Search by custom .map file
      this.filename = this.route.snapshot.queryParamMap.get("filename");
      const contourLevel = +this.route.snapshot.queryParamMap.get(
        "contourLevel"
      );
      this.biomoleculeSearchService
        .getZernikeDescriptors(emdbId, contourRepresentation)
        .then(response => {
          this.setValues(response);
        });
      this.biomoleculeSearchService
        .getSimilarBioMoleculesByMap(
          this.filename,
          contourRepresentation,
          contourLevel,
          this.volumeFilter === "On",
          minRes,
          maxRes
        )
        .then(response => {
          this.results = response.results;
          this.downloadResultFile = response.path;
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
            fill: false,
            pointRadius: 0
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
