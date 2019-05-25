const update = require("../models/updateModel");

exports.getLastUpdate = async (req, res, next) => {
  try {
    const last_update = await update.findAll({
      attributes: ["last_update"],
      limit: 1,
      order: [["last_update", "DESC"]]
    });
    if (!last_update) {
      res.status(400).send({
        message: "Last update not found"
      });
    } else {
      res.status(200).json(last_update);
    }
  } catch (err) {
    console.log(err);
    res.status(500).send({
      message: "Internal server error"
    });
  }
};
