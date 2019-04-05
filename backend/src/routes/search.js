const router = require("express").Router();
const searchController = require("../controllers/searchController");

router.get("/", (req, res, next) => {
  res.send("search");
});

router.get("/:emdbID", searchController.searchByID);

router.get(
  "/results/:emdbID/:contourShape/:isVolumeFilterOn/:minRes/:maxRes",
  searchController.searchResult
);

router.get(
  "/resultsmap/:emdbID/:contourShape/:isVolumeFilterOn/:minRes/:maxRes",
  searchController.searchResultMap
);

router.get("/zernike/:emdbID/:contourRepresentation", searchController.zernike);

router.get(
  "/zernikemap/:emdbID/:contourRepresentation",
  searchController.zernikeMap
);

module.exports = router;
