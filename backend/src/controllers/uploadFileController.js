const fs = require("fs");
const files = require("../models/filesModel");
var mkdirp = require("mkdirp");

exports.uploadFileEmMap = async (req, res, next) => {
  try {
    await files
      .build({
        file_name: req.body.filename,
        file_type: 0
      })
      .save()
      .then(fileObject => {
        saveFile(req.body.file, fileObject.file_name, fileObject.id, function(
          response
        ) {
          if (response > 0) {
            res.status(200).json(response);
          } else {
            res.status(204).json({
              msg: "ERROR: I/O error",
              details: err
            });
          }
        });
      });
  } catch (err) {
    res.status(204).json({
      msg: "ERROR: Can not upload the file",
      details: err
    });
  }
};

exports.uploadListFile = async (req, res, next) => {
  try {
    await files
      .build({
        file_name: "list.dat",
        file_type: 1
      })
      .save()
      .then(fileObject => {
        saveFile(req.body.file, fileObject.id + ".dat", fileObject.id, function(
          response
        ) {
          if (response > 0) {
            res.status(200).json(response);
          } else {
            res.status(204).json({
              msg: "ERROR: I/O error.",
              details: err
            });
          }
        });
      });
  } catch (err) {}
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
