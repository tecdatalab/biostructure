const fs = require("fs");
const files = require("../models/filesModel");
var mkdirp = require("mkdirp");

exports.uploadFileEmMap = async (req, res, next) => {
  try {
    saveFile(req.body.file, req.body.filename, "2", function(response) {
      if (response > 0) {
        res.status(200).json(response);
      } else {
        res.status(500).json({
          msg: "ERROR: I/O error",
          details: err
        });
      }
    });
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

function saveFile(file, filename, id, callback) {
  const filePath = "./files/" + id;
  mkdirp(filePath, function(err) {
    if (err) {
      callback(-1);
    } else {
      var base64Data = file.replace(
        "data:application/octet-stream;base64,",
        ""
      );
      fs.writeFile(
        filePath + "/" + filename,
        base64Data,
        { encoding: "base64" },
        function(err) {
          if (err) {
            callback(-2);
          } else {
            callback(id);
          }
        }
      );
    }
  });
}
