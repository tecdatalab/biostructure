const router = require("express").Router();
const descriptorController = require("../controllers/descriptorController");

router.get("/zernike/:emdbID", descriptorController.getZernikeDescriptors);

router.get("/zernikelist/:emdblist", descriptorController.getZernikeList);

module.exports = router;
