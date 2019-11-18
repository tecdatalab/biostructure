const update = require("../models/updateModel");

exports.getLastUpdate = async (req, res, next) => {
  try {
    const last_update = await update.findOne({
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
    res.status(500).send({
      message: "Internal server error"
    });
  }
};

exports.forceUpdate = async (req, res, next) => {
  try {
    //TO DO: Invoke updater
    console.log("update here");
    res.status(200).send({
      message: "OK"
    });
  } catch (error) {
    res.status(500).send({
      message: "Internal server error"
    });
  }
};
