import {
  async,
  ComponentFixture,
  TestBed,
  inject
} from "@angular/core/testing";
import { BiomoleculeComponent } from "./biomolecule.component";
import { BiomoleculeComparison } from "../../models/biomolecule-comparison";
import { Router } from "@angular/router";
import { Biomolecule } from "src/app/models/biomolecule";

class MockRouter {
  navigate(urls: string[], extras: string) {
    return true;
  }
}

describe("BiomoleculeComponent", () => {
  let component: BiomoleculeComponent;
  let fixture: ComponentFixture<BiomoleculeComponent>;

  beforeEach(async(() => {
    TestBed.configureTestingModule({
      declarations: [BiomoleculeComponent],
      providers: [{ provide: Router, useClass: MockRouter }]
    }).compileComponents();
  }));

  beforeEach(() => {
    fixture = TestBed.createComponent(BiomoleculeComponent);
    component = fixture.componentInstance;
    const biomoleculeComparison = new BiomoleculeComparison();
    biomoleculeComparison.biomolecule = new Biomolecule();
    biomoleculeComparison.biomolecule.id = 5555;
    biomoleculeComparison.biomolecule.full_name = "biomolecule";
    biomoleculeComparison.euc_distance = 5;
    biomoleculeComparison.biomolecule.resolution = 5.5;
    component.biomoleculeComparison = biomoleculeComparison;
    fixture.detectChanges();
  });

  it("Click on biomolecule img should call searchBiomolecule()", async(() => {
    spyOn(component, "searchBiomolecule");
    const img = fixture.debugElement.nativeElement.querySelector("#img");
    img.click();
    fixture.whenStable().then(() => {
      expect(component.searchBiomolecule).toHaveBeenCalled();
    });
  }));

  it("searchBiomolecule should call router navigate with expected params", inject(
    [Router],
    (router: Router) => {
      spyOn(router, "navigate");
      component.searchBiomolecule();
      const expectedParams = {
        queryParams: { filename: null, mapId: null, contourLevel: null },
        queryParamsHandling: "merge",
        replaceUrl: true
      };
      const expectedUrl = "result/5555";
      expect(router.navigate).toHaveBeenCalledWith(
        [expectedUrl],
        expectedParams
      );
    }
  ));
});
