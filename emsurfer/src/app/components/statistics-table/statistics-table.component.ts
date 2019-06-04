import { Component, OnInit } from "@angular/core";
import { StatisticsService } from "../../services/statistics.service";
import { Statistic } from "src/app/models/statistic";

@Component({
  selector: "app-statistics-table",
  templateUrl: "./statistics-table.component.html"
})
export class StatisticsTableComponent implements OnInit {
  statistics: Statistic[];

  constructor(private statisticsService: StatisticsService) {}

  ngOnInit() {
    this.statisticsService.getStatistics().then((statistics: Statistic[]) => {
      this.statistics = statistics;
    });
  }
}
