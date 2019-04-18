const router = require("express").Router();
const descriptorController = require("../controllers/descriptorController");

router.get("/zernike/:emdbID", descriptorController.getZernikeDescriptors);

//router.get("/zernikelist/:zernikelist", descriptorController.searchResult);

module.exports = router;
