const descriptor = require("../models/descriptorModel");
const sequelize = require("../database").sequelize;

exports.getZernikeDescriptors = async (req, res, next) => {
  try {
    const emdbid = req.params.emdbID;
    let descriptorZernike = await descriptor.findOne({
      where: {
        emd_entry_id: parseInt(emdbid)
      }
    });
    if (!descriptorZernike) {
      console.log("Biomolecule descriptors " + emdbid + " not found.");
      res.status(400).send({
        message: "Biomolecule descriptors " + emdbid + " not found."
      });
    } else {
      res.status(200).json(descriptorZernike);
    }
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};

exports.getZernikeList = async (req, res, next) => {};
