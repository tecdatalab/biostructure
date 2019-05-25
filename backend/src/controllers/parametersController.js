const parameters = require("../parameters.json");

exports.getParameters = async (req, res, next) => {
  res.status(200).json(parameters);
};
