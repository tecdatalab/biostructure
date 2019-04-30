const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const benchmark = sequelize.define(
  "benchmark_history",
  {
    id: {
      type: Sequelize.INTEGER,
      primaryKey: true
    },
    date_time: {
      type: Sequelize.DATE
    },
    ip: {
      type: Sequelize.INTEGER
    },
    user_id: {
      type: Sequelize.INTEGER
    },
    representation_id: {
      type: Sequelize.INTEGER
    },
    volume_filter_id: {
      type: Sequelize.INTEGER
    },
    top_results: {
      type: Sequelize.INTEGER
    },
    emd_list: {
      type: Sequelize.JSON
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "benchmark_history"
  }
);

module.exports = benchmark;
