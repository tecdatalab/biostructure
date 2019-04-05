const router = require("express").Router();
const searchController = require("../controllers/searchController");

router.get("/", (req, res, next) => {
  res.send("search");
});

router.get(
  '/:emdbID',
  searchController.searchByID
);

router.get(
  '/results/:emdbID/:isVolumeFilterOn/:minRes/:maxRes',
  searchController.searchResult
);

router.get(
  '/zernike/:emdbID/:contourRepresentation',
  searchController.zernike
);

module.exports = router;
