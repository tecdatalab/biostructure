const search_history = require("../models/searchHistoryModel");
const Op = require("../database").Op;

exports.getSearchHistory = (req, res, next) => {
    const reqUser = req.authenticatedUser;
    try {
        const emdbid = req.params.userID;
        let searches = await search_history.findAll({
            where: {
                id: reqUser.id
            }
        })
    } catch (err) {
        res.status(500).send({
            message: "Internal server error"
        });
    }
}