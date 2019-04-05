import { Injectable } from '@angular/core';

@Injectable({
  providedIn: 'root'
})
export class ContourRepresentationService {
  constructor() {}

  getContourShapes() {
    const shapes = [
      'EMDB contour',
      'EMDB contour + 1/3 core',
      'EMDB contour + 2/3 core',
      'EMDB contour + 1/3 + 2/3 core',
      'EMDB contour + 1 std dev'
    ];
    return shapes;
  }
}
