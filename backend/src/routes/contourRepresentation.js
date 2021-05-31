const router = require("express").Router();
const contourController = require("../controllers/contourRepresentationController");

router.get(
  "/",
  contourController.getContourRepresentation
);

module.exports = router;
