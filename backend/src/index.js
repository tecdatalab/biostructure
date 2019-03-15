const express = require('express')
const app = express()
const indexRoutes = require('./routes/index');

//settings 
app.set('port', process.env.PORT || 3000)

//middlewares
app.use(express.json());
app.use(express.urlencoded({extended: false}));

//routes
app.use(indexRoutes)

app.listen(app.get('port'), () => {
    console.log('Server on port ', app.get('port'))
});