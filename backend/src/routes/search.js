const router = require("express").Router();
const searchController = require("../controllers/searchController");
const userController = require("../controllers/userController");

router.get("/:app/:ID", searchController.searchByID);

router.get(
  "/results/:emdbID/:contourRepresentation/:isVolumeFilterOn/:minRes/:maxRes/:typeDescriptor/:topCan",
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

router.get("/:ID",searchController.getCathDetail);

router.get("/getResultsPDB/:ID/:rep/:type/:cath/:len/:top", searchController.getResultsPDB)

module.exports = router;
