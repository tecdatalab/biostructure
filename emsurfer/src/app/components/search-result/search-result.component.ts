import { Component, OnInit, ViewChild, ElementRef } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { Chart } from 'chart.js';
import { BiomoleculeSearchService } from '../../services/biomolecule-search.service';
import { Biomolecule } from 'src/app/models/biomolecule';
import { BiomoleculeComparison } from 'src/app/models/biomolecule-comparison';

@Component({
  selector: 'app-search-result',
  templateUrl: './search-result.component.html',
  providers: [BiomoleculeSearchService]
})
export class SearchResultComponent implements OnInit {
  @ViewChild('canvas') canvasElementRef: ElementRef;
  chart: Chart;
  biomolecule: Biomolecule;
  filename: string;
  results: BiomoleculeComparison[];
  volumeFilter: string;
  isSearchById: boolean;
  descriptors = [];
  values = [];

  constructor(
    private biomoleculeSearchService: BiomoleculeSearchService,
    private route: ActivatedRoute
  ) {}

  ngOnInit() {
    this.biomolecule = new Biomolecule();
    this.route.params.subscribe(params => {
      this.biomolecule.emdb_id = params.emdbId;
      this.load();
    });
  }

  private load() {
    const emdbId = +this.route.snapshot.paramMap.get('emdbId') || null;
    const contourRepresentation = +this.route.snapshot.queryParamMap.get(
      'contourRepresentation'
    );
    const minRes = +this.route.snapshot.queryParamMap.get('minRes');
    const maxRes = +this.route.snapshot.queryParamMap.get('maxRes');
    this.volumeFilter = this.route.snapshot.queryParamMap.get('volumeFilter');
    if (emdbId) {
      this.biomolecule = this.biomoleculeSearchService.getBiomolecule(emdbId);
      this.values = this.biomoleculeSearchService.getZernikeDescriptors(
        emdbId,
        contourRepresentation
      );
      this.isSearchById = true;
    } else {
      this.filename = this.route.snapshot.queryParamMap.get('filename');
      this.values = this.biomoleculeSearchService.getZernikeDescriptors(
        emdbId,
        contourRepresentation
      );
    }
    this.descriptors = Array.from(
      new Array(this.values.length),
      (val, index) => index + 1
    ); // [1,2,3...N]
    this.results = this.biomoleculeSearchService.getSimilarBioMolecules(
      5555,
      contourRepresentation,
      this.volumeFilter === 'On',
      minRes,
      maxRes
    );
    const context = this.canvasElementRef.nativeElement;
    this.initChart(context);
  }

  private initChart(context: ElementRef) {
    this.chart = new Chart(context, {
      type: 'line',
      data: {
        labels: this.descriptors,
        datasets: [
          {
            data: this.values,
            borderColor: 'black',
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
                labelString: 'Zernike Descriptor Number',
                fontSize: 24
              }
            }
          ],
          yAxes: [
            {
              display: true,
              scaleLabel: {
                display: true,
                labelString: 'Value',
                fontSize: 24
              }
            }
          ]
        }
      }
    });
  }
}
