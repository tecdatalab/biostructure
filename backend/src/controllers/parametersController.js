fs = require("fs");

exports.getParameters = async (req, res, next) => {
  fs.readFile("./src/parameters.json", "utf8", (err, parameters) => {
    if (err) {
      res.status(500).send({
        message: "Internal server error"
      });
    }
    res.status(200).json(JSON.parse(parameters));
  });
};

exports.setParameters = async (req, res, next) => {
  try {
    var data = {
      volume_filter_min: isNaN(req.params.volumeMin)
        ? 0
        : parseFloat(req.params.volumeMin),
      volume_filter_max: isNaN(req.params.volumeMax)
        ? 0
        : parseFloat(req.params.volumeMax),
      hits_number: isNaN(req.params.hits) ? 10 : parseFloat(req.params.hits),
      update_rate: isNaN(req.params.update)
        ? 100
        : parseFloat(req.params.update)
    };
    fs.writeFile("./src/parameters.json", JSON.stringify(data), err => {
      if (err) {
        res.status(500).send({
          message: "Internal server error"
        });
      } else {
        res.status(200).send({
          message: "OK"
        });
      }
    });
  } catch (error) {
    res.status(500).send({
      message: "Internal server error"
    });
  }
};
