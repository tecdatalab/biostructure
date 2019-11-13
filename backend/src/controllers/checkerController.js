const biomolecule = require("../models/biomoleculeModelEMD");

exports.checkBiomolecule = async (req, res, next) => {
  emdbid = req.params.emdbID;
  try {
    let biomolecules = await biomolecule.count({
      where: {
        id: parseInt(emdbid)
      }
    });
    if (!biomolecules) {
      console.log("Biomolecule " + emdbid + " not found.");
      res.status(400).send({
        message: "Biomolecule " + emdbid + " not found."
      });
    } else {
      res.status(200).json(biomolecules);
    }
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};
