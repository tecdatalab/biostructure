const router = require('express').Router();

router.get('/', (req, res, next) => {
    res.send('search');
});

module.exports = router