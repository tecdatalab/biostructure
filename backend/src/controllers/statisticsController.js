const biomolecule = require("../models/biomoleculeModelEMD");
const Sequelize = require("sequelize");
const Op = require("../database").Op;

exports.getStatistics = async (req, res) => {
  try {
    const allRes = await biomolecule.count();
    const reslet5 = await biomolecule.count({
      where: { resolution: { [Op.lte]: 5 } }
    });
    const resgt5let10 = await biomolecule.count({
      where: { resolution: { [Op.and]: [{ [Op.gt]: 5 }, { [Op.lte]: 10 }] } }
    });
    const resgt10let15 = await biomolecule.count({
      where: { resolution: { [Op.and]: [{ [Op.gt]: 10 }, { [Op.lte]: 15 }] } }
    });
    const resgt15let20 = await biomolecule.count({
      where: { resolution: { [Op.and]: [{ [Op.gt]: 15 }, { [Op.lte]: 20 }] } }
    });
    const resgt20 = await biomolecule.count({
      where: { resolution: { [Op.gt]: 20 } }
    });
    const statistics = [
      { name: "All Entries with resolution data", value: allRes },
      { name: "Entries with resolution <=5Å", value: reslet5 },
      { name: "Entries with resolution >5Å and <=10Å", value: resgt5let10 },
      { name: "Entries with resolution >10Å and <=15Å", value: resgt10let15 },
      { name: "Entries with resolution >15Å and <=20Å", value: resgt15let20 },
      { name: "Entries with resolution >20Å", value: resgt20 }
    ];
    res.status(200).json(statistics);
  } catch (err) {
    res.status(500).send({
      message: "Backend error"
    });
  }
};
