const router = require("express").Router();
const descriptorController = require("../controllers/descriptorController");

router.get(
  "/zernike/:emdbID/:contourRepresentation",
  descriptorController.getZernikeDescriptors
);

router.get(
  "/zernikelist/:emdblist/:contourRepresentation",
  descriptorController.getZernikeList
);

module.exports = router;
