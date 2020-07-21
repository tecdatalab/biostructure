const router = require("express").Router();
const visualizerController = require("../controllers/visualizerController");

router.get(
  "/:emdId",
  visualizerController.searchByID
);

module.exports = router;