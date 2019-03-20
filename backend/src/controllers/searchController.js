const biomolecule = require('../models/biomolecule')

exports.searchByID = async (req, res, next) => {
    emdbid = req.params.emdbID
    try{
        let biomolecules = await biomolecule.findOne({
            where: {
                id: parseInt(emdbid)
            }
        })
        res.status(200).json(biomolecules);
    } catch (error) {
        console.log(error)
        res.status(500).json({
            msg: 'Biomolecule not found',
            details: err
        })
    }
  }