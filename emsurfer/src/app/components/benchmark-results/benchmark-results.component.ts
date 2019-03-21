import { Component, OnInit } from '@angular/core';
import { ActivatedRoute } from '@angular/router';
import { BiomoleculeSearchService } from 'src/app/services/biomolecule-search.service';
import { CustomFile } from 'src/app/models/custom-file';
import { FileDownloadService } from 'src/app/services/file-download.service';

@Component({
  selector: 'app-benchmark-results',
  templateUrl: './benchmark-results.component.html'
})
export class BenchmarkResultsComponent implements OnInit {
  files: CustomFile[];
  compressedFilePath: string;
  constructor(
    private route: ActivatedRoute,
    private biomoleculeSearchService: BiomoleculeSearchService,
    private fileDownloadService: FileDownloadService
  ) {}

  ngOnInit() {
    const contourRepresentationId = +this.route.snapshot.queryParamMap.get(
      'contourRepresentation'
    );
    const volumeFilter =
      this.route.snapshot.queryParamMap.get('volumeFilter') === 'On';
    const topResults = +this.route.snapshot.queryParamMap.get('topResults');
    const emdbIdList = +this.route.snapshot.queryParamMap.get('emdbIdList');
    const emdbIdListFile = +this.route.snapshot.queryParamMap.get('emdbIdListFile');
    if (emdbIdList){
      this.files = this.biomoleculeSearchService.getBatchBiomolecules(
        emdbIdList,
        contourRepresentationId,
        volumeFilter,
        topResults
      );
    } else if(emdbIdListFile) {
      this.files = this.biomoleculeSearchService.getBatchBiomoleculesByFileId(
        emdbIdListFile,
        contourRepresentationId,
        volumeFilter,
        topResults
      );
    }
    // don't forget to include the id of the compressed file created in the server
    this.compressedFilePath = this.fileDownloadService.getBenchmarkResultCompressedFilePath(1444);
  }
}
