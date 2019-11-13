const contourRepresentation = require("../models/contourRepresentationModel");

exports.getContourRepresentation = async (req, res, next) => {
  try {
    let contours = await contourRepresentation.findAll({});
    if (!contours) {
      console.log("Contours Representation not found.");
      res.status(400).send({
        message: "Contours Representation not found."
      });
    } else {
      res.status(200).json(contours);
    }
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};
