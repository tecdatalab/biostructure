const router = require("express").Router();
const searchController = require("../controllers/searchController");

router.get("/", (req, res, next) => {
  res.send("search");
});

router.get(
  '/:emdbID/:isVolumeFilterOn/:minRes/:maxRes',
  searchController.searchByID
);

module.exports = router;
