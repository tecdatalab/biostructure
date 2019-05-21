const router = require("express").Router();
const searchController = require("../controllers/searchController");
const userController = require("../controllers/userController");

router.get("/:emdbID", searchController.searchByID);

router.get(
  "/results/:emdbID/:contourRepresentation/:isVolumeFilterOn/:minRes/:maxRes",
  userController.isUserLoggedIn,
  searchController.saveSearch,
  searchController.searchResult
);

router.get(
  "/resultsmap/:filename/:contourRepresentation/:contourLevel/:isVolumeFilterOn/:minRes/:maxRes",
  userController.isUserLoggedIn,
  searchController.saveSearch,
  searchController.searchResultMap
);

router.get("/zernike/:emdbID/:contourRepresentation", searchController.zernike);

router.get(
  "/zernikemap/:emdbID/:contourRepresentation",
  searchController.zernikeMap
);

module.exports = router;
