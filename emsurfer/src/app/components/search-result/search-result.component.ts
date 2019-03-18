import { Component, OnInit, ViewChild, ElementRef } from '@angular/core';
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
  results: BiomoleculeComparison[];
  descriptors = [];
  values = [];

  constructor(private biomoleculeSearchService: BiomoleculeSearchService) { }

  ngOnInit() {
    this.biomolecule = this.biomoleculeSearchService.getBiomolecule(1413);
    this.values = this.biomoleculeSearchService.getZernikeDescriptors(1413);
    this.descriptors = Array.from(
      new Array(this.values.length),
      (val, index) => index + 1
    ); // [1,2,3...N]
    const context = this.canvasElementRef.nativeElement;
    this.initChart(context);
    this.results = this.biomoleculeSearchService.getSimilarBioMolecules(1413);
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
