import { Component, OnInit } from "@angular/core";
import { SearchHistory } from "src/app/models/search-history";
import { SearchHistoryService } from "src/app/services/search-history.service";

@Component({
  selector: "app-search-history",
  templateUrl: "./search-history.component.html",
  styleUrls: ["./search-history.component.css"]
})
export class SearchHistoryComponent implements OnInit {
  constructor(private searchHistoryService: SearchHistoryService) {}
  records: SearchHistory[];

  ngOnInit() {
    this.searchHistoryService
      .getSearchHistory()
      .then((response: SearchHistory[]) => {
        console.log(response);
        this.records = response;
      });
  }
}
