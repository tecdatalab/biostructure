const biomolecule = require('../models/biomolecule')

exports.searchByID = async (req, res, next) => {
    emdbid = req.params.emdbID
    try{
        let biomolecules = await biomolecule.findAll({
            attributes: ['id'],
            where: {
                id: parseInt(emdbid)
            }
        })
        res.json({
            result: 'OK',
            data: biomolecules,
            message: 'Biomolecule found'
        });
    } catch (error) {
        console.log(error)
        res.json({
            result: 'Failed',
            data: [],
            message: 'Biomolecule not found'
        })
    }
  }