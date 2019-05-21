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
  minDate = new Date();
  currentPage = -1;
  rowsPerPage = 6;

  searchByEmdbID(record) {
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

  nextPage() {
    this.currentPage += this.rowsPerPage - 1;
  }

  previousPage() {
    this.currentPage -= this.rowsPerPage - 1;
  }

  filterFunction(records) {
    return records.filter(record => {
      //console.log(new Date(record["date_time"]));
      if (new Date(record["date_time"]) < this.minDate) {
        return true;
      }
      return false;
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
