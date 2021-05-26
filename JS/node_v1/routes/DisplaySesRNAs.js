
// Initializing express server 
const express = require('express')
// const router = express.Router()
const router = express()
// Defining output port 
const port = 9001

// Spawn child process 
const {spawn} = require('child_process')

router.get('/', function(req, res, next) {

 var dataToSend;
 // Spawn new child process to call the python script
 // Returns information on splice variant given species and gene names 
 const pyFunction_spliceVariant = '/home/user1/Dropbox/Research/Neurobiology_PhD/Huang/Projects/CellReadR_Temp/Code/kCellReadR/json_spliceInformation.py'
 const python = spawn('python', [pyFunction_spliceVariant, 'Rat', 'Bcl11b'])
//  const python = spawn('python', ['test.py'])

 // Collect data from script
 python.stdout.on('data', function (data) {
  console.log('Pipe data from python script ...')
  console.log('pyFunction_call')
  console.log(`Data:${data}`);
  dataToSend = data.toString();
//   dataToSend.push(data)
//   dataToSend = data
 })

 // Close child python process 
 python.on('close', (code) => {
 console.log(`child process close all stdio with code ${code}`)
 // Send data to browser
 res.send(dataToSend)
 })
})


// Listen 
router.listen(port, () => console.log(`DisplaySesRNAs app listening on port ${port}!`))

module.exports = router
