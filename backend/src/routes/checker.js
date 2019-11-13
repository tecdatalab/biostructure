const router = require("express").Router();
const checkerController = require("../controllers/checkerController");

router.get(
  "/emdbID/:emdbID",
  checkerController.checkBiomolecule
);

module.exports = router;
