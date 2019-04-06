const contourRepresentation = require("../models/contourRepresentationModel");
const sequelize = require("../database").sequelize;

exports.getContourRepresentation = async (req, res, next) => {
  try {
    let contours = await contourRepresentation.findAll({});
    if (!contours) {
      console.log("Contours Representation not found.");
      res.status(500).json({
        msg: "Contours Representation not found.",
        details: err
      });
    } else {
      res.status(200).json(contours);
    }
  } catch (err) {
    res.status(204).json({
      msg: "Backend error",
      details: err
    });
  }
};
