import { Component, OnInit } from "@angular/core";
import { Router } from "@angular/router";
import { SearchHistory } from "src/app/models/search-history";
import { UserService } from "src/app/services/user.service";
import { SearchHistoryService } from "src/app/services/search-history.service";

@Component({
  selector: "app-search-history",
  templateUrl: "./search-history.component.html",
  styleUrls: ["./search-history.component.css"]
})
export class SearchHistoryComponent implements OnInit {
  constructor(
    private userService: UserService,
    private searchHistoryService: SearchHistoryService,
    private router: Router
  ) {}
  records: SearchHistory[];
  minDate: Date;
  maxDate: Date;
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
      if (
        (this.minDate == null ||
          new Date(this.minDate) <= new Date(record["date_time"])) &&
        (this.maxDate == null ||
          new Date(record["date_time"]) <= new Date(this.maxDate))
      ) {
        return true;
      }
      return false;
    });
  }

  ngOnInit() {
    if (this.userService.isUserLoggedIn()) {
      this.searchHistoryService
        .getSearchHistory()
        .then((response: SearchHistory[]) => {
          this.records = response;
        });
    } else {
      this.router.navigate(["/home"]);
    }
  }
}
