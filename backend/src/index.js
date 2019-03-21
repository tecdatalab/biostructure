const express = require('express');
const cors = require('cors');
const app = express();
const bodyParser = require('body-parser')
const expressValidator = require('express-validator');
const indexRoutes = require('./routes/index');
const searchRoutes = require('./routes/search');
const uploadFileRoutes = require("./routes/uploadFile");

//settings 
app.set('port', process.env.PORT || 3001)

//middlewares
app.use(bodyParser.urlencoded({limit: '50mb', extended: true }))
app.use(bodyParser.json({ limit: '50mb', extended: true }))
app.use(express.json());
app.use(cors())
//app.use(express.urlencoded({extended: false}));
app.use(expressValidator())

//routes
app.use(indexRoutes);
app.use('/img', express.static('public/img'));
app.use('/results', express.static('public/results'));
app.use('/search', searchRoutes);
app.use('/upload', uploadFileRoutes);

app.listen(app.get('port'), () => {
    console.log('Server on port ', app.get('port'))
});