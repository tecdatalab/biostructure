const search_history = require("../models/searchHistoryModel");
const Op = require("../database").Op;

exports.getSearchHistory = async (req, res, next) => {
  const reqUser = req.authenticatedUser;
  try {
    let searches = await search_history.findAll({
      where: {
        user_id: reqUser.id
      }
    });
    res.status(200).json(searches);
  } catch (err) {
    console.log(err);
    res.status(500).send({
      message: "Internal server error"
    });
  }
};
