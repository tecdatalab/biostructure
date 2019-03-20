import { Injectable } from '@angular/core';

@Injectable({
  providedIn: 'root'
})
export class FileUploadService {

  constructor() { }

  uploadIdListFile(base64File: string) {
    return 5555;
  }

  uploadIdListString(idList: string){
    return 1234;
  }

  uploadEmMap(base64File: string){
    return 4444;
  }
}
