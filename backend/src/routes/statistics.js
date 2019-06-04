const router = require("express").Router();
const statisticsController = require("../controllers/statisticsController");

router.get("/", statisticsController.getStatistics);

module.exports = router;
