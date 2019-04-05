import { Component, OnInit } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { BiomoleculeSearchService } from 'src/app/services/biomolecule-search.service';
import { CustomFile } from 'src/app/models/custom-file';

@Component({
  selector: 'app-benchmark-results',
  templateUrl: './benchmark-results.component.html'
})
export class BenchmarkResultsComponent implements OnInit {
  files: CustomFile[];
  constructor(
    private route: ActivatedRoute,
    private biomoleculeSearchService: BiomoleculeSearchService
  ) {}

  ngOnInit() {
    const contourRepresentationId = +this.route.snapshot.queryParamMap.get(
      'contourRepresentation'
    );
    const volumeFilter =
      this.route.snapshot.queryParamMap.get('volumeFilter') === 'On';
    const topResults = +this.route.snapshot.queryParamMap.get('topResults');
    this.files = this.biomoleculeSearchService.getBatchBiomolecules(
      1,
      contourRepresentationId,
      volumeFilter,
      topResults
    );
  }
}
