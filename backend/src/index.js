const express = require('express')
const app = express()
const expressValidator = require('express-validator')
const indexRoutes = require('./routes/index');
const searchRoutes = require('./routes/search')

//settings 
app.set('port', process.env.PORT || 3000)

//middlewares
app.use(express.json());
app.use(express.urlencoded({extended: false}));
app.use(expressValidator())

//routes
app.use(indexRoutes);
app.use('/img', express.static('public/img'));
app.use('/search', searchRoutes);

app.listen(app.get('port'), () => {
    console.log('Server on port ', app.get('port'))
});