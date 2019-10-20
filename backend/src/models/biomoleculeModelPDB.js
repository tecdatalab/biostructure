const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;
const cathInfo = require("cathModel");

const biomolecule = sequelize.define(
    "pdb_entry",
    {
        id: {
            type: Sequelize.INTEGER,
            primaryKey : true
        },
        pdb_id: {
            type: Sequelize.STRING
        },
        numbers_descriptor: {
            type: Sequelize.JSON
        },
        secuence: {
            type: Sequelize.STRING
        },
        png_img_3d: {
            type: Sequelize.STRING
        },
        gif_img_3d: {
            type: Sequelize.STRING
        }
    },
    {
        timestamps: false,
        freezeTableName: true,
        tableName: "atomic_structure"
    }
);

//biomolecule.findAll({
//    include: [{
//        model: Tool,
//        required: true,
//        right: true // has no effect, will create an inner join
//    }]
//});
module.exports = biomolecule;
