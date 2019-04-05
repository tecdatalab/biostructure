const router = require('express').Router();

router.get('/', (req, res, next) => {
    res.send('connection to the server');
});

module.exports = router