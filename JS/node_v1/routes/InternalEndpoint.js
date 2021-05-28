// Initializing express server 
const express = require('express')
// const router = express.Router()
const router = express()
// Defining output port 
const port = 9000

// Need to first parse our incoming requests
router.use(express.json());
router.use(express.urlencoded({ extended: true }));

// Spawn child process 
const {spawn} = require('child_process')

// Fetches splice variant information 
// Given species and gene name
router.post('/POST_spliceVariants', function(req, res, next) {
  console.log(`Test 1:${Object.keys(req)}`);
  console.log(`Species:${req.body.species}`);
  console.log(`Gene:${req.body.gene}`);

  // Spawn new child process to call the python script
  // Returns information on splice variant given species and gene names 
  const pyFunction_spliceVariant = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR_Temp/Code/kCellReadR/json_spliceInformation.py'
  const parameters = [pyFunction_spliceVariant, req.body.species, req.body.gene]
  const python = spawn('python', parameters)

  // Declaring global variable prior to python function 
  var dataToSend;
  // Collect data from script
  python.stdout.on('data', function (data) {
  console.log('Pipe data from python script ...')
  // Global scope
  dataToSend = data.toString();
  })

  // Close child python process 
  python.on('close', (code) => {
  console.log(`child process close all stdio with code ${code}`)
  // Send data to browser
  res.send(dataToSend)
  })
})

// Fetches candidate sesRNAs 
// Given species, gene name, and sesRNA selection parameters 
router.post('/POST_sesRNAs', function(req, res, next) {
  console.log(`Test 1:${Object.keys(req)}`);
  console.log(`Test:${req.body.species}`);
  console.log(`Test:${req.body.gene}`);
  console.log(`Test:${req.body.minGC}`);
  console.log(`Test:${req.body.choice_ATG}`);

  // Defining function name/path and parameters for generating sesRNAs
  const pyFunction_sesRNAs = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR_Temp/Code/kCellReadR/json_sesRNAs.py'
  const parameters = [pyFunction_sesRNAs, 
                      req.body.species, 
                      req.body.gene,
                      req.body.spliceVariant,
                      req.body.variantTable, 
                      req.body.searchSeq,
                      req.body.seqDirection,
                      req.body.len_sesRNA,
                      req.body.minTGG,
                      req.body.maxStop, 
                      req.body.choice_ATG,
                      req.body.minGC, 
                      req.body.maxGC, 
                      req.body.dist_cTGG,
                      req.body.dist_stop_cTGG]


  // Spawn new child process to call the python script
  // Generates sesRNAs
  const python = spawn('python', parameters) 

  // Declaring global variable prior to python function 
  var dataToSend;
  // Collect data from script
  python.stdout.on('data', function (data) {
  console.log('Pipe data from python script for sesRNAs ...')
  // Global scope
  dataToSend = data.toString();
  console.log(dataToSend)
  })

  // Close child python process 
  python.on('close', (code) => {
  console.log(`child process close all stdio with code ${code}`)
  // Send data to browser
  console.log(dataToSend)
  res.send(dataToSend)
  })
})

// Listen 
router.listen(port, () => console.log(`DisplaySesRNAs app listening on port ${port}!`))

module.exports = router