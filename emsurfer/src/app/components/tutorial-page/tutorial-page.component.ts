import { Component, OnInit } from '@angular/core';
import { interval as observableInterval } from "rxjs";
import { takeWhile, scan, tap } from "rxjs/operators";

@Component({
  selector: 'app-tutorial-page',
  templateUrl: './tutorial-page.component.html',
  styleUrls: ['./tutorial-page.component.css'],
})
export class TutorialPageComponent implements OnInit {

  constructor() { }

  ngOnInit() {
  }

  ScrollTop(){
    window.scroll(0,300);
  }

}
