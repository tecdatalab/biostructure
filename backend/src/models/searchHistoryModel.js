const Sequelize = require("sequelize");
const sequelize = require("../database").sequelize;

const search_history = sequelize.define(
  "search_history",
  {
    id: {
      type: Sequelize.INTEGER,
      primaryKey: true
    },
    date_time: {
      type: Sequelize.DATE
    },
    ip: {
      type: Sequelize.STRING
    },
    emd_entry_id: {
        type: Sequelize.INTEGER
    },
    name_file: {
        type: Sequelize.STRING
    },
    counter_level: {
        type: Sequelize.DOUBLE
    },
    representation_id: {
        type: Sequelize.INTEGER
    },
    volume_filter_id: {
        type: Sequelize.INTEGER
    },
    resolution_filter_min: {
        type: Sequelize.DOUBLE
    },
    resolution_filter_max: {
        type: Sequelize.DOUBLE
    }
  },
  {
    timestamps: false,
    freezeTableName: true,
    tableName: "search_history"
  }
);

module.exports = search_history;