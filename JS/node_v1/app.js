var createError = require('http-errors');
var express = require('express');
var path = require('path');
var cookieParser = require('cookie-parser');
var logger = require('morgan');

// For cors 
const cors = require("cors")

// Requiring express to use user defined API/route 
var router_sesRNAs = require("./routes/DisplaySesRNAs")
var router_spliceVariantInfo = require("./routes/SpliceVariantInfo")

var app = express();

// view engine setup
app.set('views', path.join(__dirname, 'views'));
app.set('view engine', 'jade');

app.use(logger('dev'));
// app.use(express.urlencoded({ extended: false }));
// app.use(express.urlencoded())
app.use(express.urlencoded({ extended: true }));
app.use(express.json()) // for handling any json requests 

app.use(cookieParser());
app.use(express.static(path.join(__dirname, 'public')));
app.use(express.json()); //Used to parse JSON bodies
app.use(express.urlencoded()); //Parse URL-encoded bodies

// For cors 
app.use(cors())

// Telling express to use user definied route ... 
app.use("/DisplaySesRNAs", router_sesRNAs)
app.use("/DisplaySesRNAs", router_spliceVariantInfo)

// catch 404 and forward to error handler
app.use(function(req, res, next) {
  next(createError(404));
});

// error handler
app.use(function(err, req, res, next) {
  // set locals, only providing error in development
  res.locals.message = err.message;
  res.locals.error = req.app.get('env') === 'development' ? err : {};

  // render the error page
  res.status(err.status || 500);
  res.render('error');
});

module.exports = app;
