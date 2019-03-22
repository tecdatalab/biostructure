const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;
const Op = require("../database").Op;

const benchmark = sequelize.define(
  "benchmark",
  {
    id: {
      type: Sequelize.INTEGER,
      primaryKey: true
    },
    date: {
      type: Sequelize.DATE
    },
    fk_user: {
      type: Sequelize.INTEGER
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "benchmark"
  }
);

module.exports = benchmark;
