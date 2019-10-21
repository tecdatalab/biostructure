const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const cath = sequelize.define(
    "cathDetail",
    {
        atomic_structure_id: {
            type: Sequelize.INTEGER,
            primaryKey : true
        },
        class_number: {
            type: Sequelize.INTEGER
        },
        architecture_number: {
            type: Sequelize.INTEGER
        },
        topology_number: {
            type: Sequelize.INTEGER
        },
        homologous_superfamily_number: {
            type: Sequelize.INTEGER
        },
        domain_length: {
            type: Sequelize.INTEGER
        }
    },
    {
        timestamps: false,
        freezeTableName: true,
        tableName: "cath_atomic_structure"
    }
);

module.exports = cath;