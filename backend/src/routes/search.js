const router = require("express").Router();
const searchController = require("../controllers/searchController");

router.get("/:emdbID", searchController.searchByID);

router.get(
  "/results/:emdbID/:contourRepresentation/:isVolumeFilterOn/:minRes/:maxRes",
  searchController.searchResult
);

router.get(
  "/resultsmap/:filename/:contourRepresentation/:contourLevel/:isVolumeFilterOn/:minRes/:maxRes",
  searchController.searchResultMap
);

router.get("/zernike/:emdbID/:contourRepresentation", searchController.zernike);

router.get(
  "/zernikemap/:emdbID/:contourRepresentation",
  searchController.zernikeMap
);

module.exports = router;
