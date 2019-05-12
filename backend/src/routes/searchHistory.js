const searchHistoryController = require("../controllers/searchHistoryController");
const userController = require("../controllers/userController");
const router = require("express").Router();

router.get(
  "/",
  userController.verifyUserToken,
  searchHistoryController.getSearchHistory
);
