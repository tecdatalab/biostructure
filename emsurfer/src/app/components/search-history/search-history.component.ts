import { Component, OnInit } from "@angular/core";
import { Router } from "@angular/router";
import { SearchHistory } from "src/app/models/search-history";
import { SearchHistoryService } from "src/app/services/search-history.service";

@Component({
  selector: "app-search-history",
  templateUrl: "./search-history.component.html",
  styleUrls: ["./search-history.component.css"]
})
export class SearchHistoryComponent implements OnInit {
  constructor(
    private searchHistoryService: SearchHistoryService,
    private router: Router
  ) {}
  records: SearchHistory[];

  searchByEmdbID(record) {
    console.log(record);
    const url = "result/" + record.emd_entry_id;
    const params = {
      contourRepresentation: record.representation_id,
      volumeFilter: record.volume_filter_id,
      minResolution: record.resolution_filter_min,
      maxResolution: record.resolution_filter_max
    };
    this.router.navigate([url], {
      queryParams: params
    });
  }

  ngOnInit() {
    this.searchHistoryService
      .getSearchHistory()
      .then((response: SearchHistory[]) => {
        console.log(response);
        this.records = response;
      });
  }
}
