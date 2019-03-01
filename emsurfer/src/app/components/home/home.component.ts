import { ActivatedRoute } from '@angular/router'
import { Component, OnInit } from '@angular/core';

@Component({
  selector: 'app-home',
  templateUrl: './home.component.html'
})
export class HomeComponent implements OnInit {

  constructor(private route: ActivatedRoute) { }

  params = {
    id: '',
    volume: '',
    min: '',
    max: ''
  }

  ngOnInit() {
    this.route.paramMap
      .subscribe(params => {
       this.params.id  = params.get('id');
       this.params.volume  = params.get('volume');
       this.params.min = params.get('min');
       this.params.max = params.get('max');
      });
  }
}
