const Sequelize = require('sequelize');
const sequelize = require('../database').sequelize;
const Op = require('../database').Op;

const biomolecule = sequelize.define('emd_entry', {
  id: {
      type: Sequelize.INTEGER,
      primaryKey: true
  },
  full_name: {
      type: Sequelize.STRING
  },
  acronym: {
      type: Sequelize.STRING
  },
  volume: {
      type: Sequelize.DOUBLE
  },
  resolution: {
      type: Sequelize.DOUBLE
  },
  image_url: {
      type: Sequelize.STRING
  },
  xml_url: {
      type: Sequelize.STRING
  },
  map_url: {
      type: Sequelize.STRING
  }
},{
  timestamps: false,
  freezeTableName: true,
  tableName: 'emd_entry'
});

module.exports = biomolecule