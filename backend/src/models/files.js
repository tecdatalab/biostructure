const Sequelize = require('sequelize');
const sequelize = require('../database').sequelize;
const Op = require('../database').Op;

const files = sequelize.define('files', {
  id: {
      type: Sequelize.INTEGER,
      primaryKey: true
  },
  file_name: {
      type: Sequelize.STRING
  },
  file_type: {
      type: Sequelize.INTEGER
  }
},{
  timestamps: false,
  freezeTableName: true,
  tableName: 'files'
});

module.exports = files