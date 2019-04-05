const fs = require("fs");
const files = require("../models/files");
const sequelize = require("../database").sequelize;
var mkdirp = require("mkdirp");

exports.uploadFileEmMap = async (req, res, next) => {
  let fileId;
  try {
    await files
      .build({
        file_name: req.body.filename,
        file_type: 0
      })
      .save()
      .then(fileObject => {
        const filePath = "./files/" + fileObject.dataValues.id;
        fileId = fileObject.dataValues.id;
        mkdirp(filePath, function(err) {
          if (err) {
            res.status(204).json({
              msg: "Path already exits.",
              details: err
            });
          } else {
            console.log(fileObject.dataValues.id);
            var base64Data = req.body.file.replace(
              "data:application/octet-stream;base64,",
              ""
            );
            fs.writeFile(
              filePath + "/" + req.body.filename,
              base64Data,
              { encoding: "base64" },
              function(err) {
                console.log("File created");
              }
            );
          }
        });
      });
    return res.status(200).json(fileId);
  } catch (err) {
    res.status(204).json({
      msg: "Can not upload the file",
      details: err
    });
  }
};
