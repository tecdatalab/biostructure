const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const contourRepresentation = sequelize.define(
  "contourRepresentation",
  {
    id: {
      type: Sequelize.INTEGER,
      primaryKey: true
    },
    name: {
      type: Sequelize.STRING
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "representation"
  }
);

module.exports = contourRepresentation;
