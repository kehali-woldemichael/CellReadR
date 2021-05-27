import React, { Component } from 'react'
import {MuiThemeProvider} from "@material-ui/core/styles"

import { TextField } from '@material-ui/core'

import { Button } from '@material-ui/core'

import AppBar from '@material-ui/core/AppBar'
import { Toolbar } from '@material-ui/core'
import { Typography } from '@material-ui/core'
 
// import { MenuIcon } from '@material-ui/icons/Menu'
// import { InputLabel } from '@material-ui/core'
// import { IconButton } from '@material-ui/core'

import Select from '@material-ui/core/Select'
import InputLabel from '@material-ui/core/InputLabel'
import Grid from '@material-ui/core/Grid'
import Box from '@material-ui/core/Box';
import MenuItem from '@material-ui/core/MenuItem';


import Table from '@material-ui/core/Table';
import TableBody from '@material-ui/core/TableBody';
import TableCell from '@material-ui/core/TableCell';
import TableContainer from '@material-ui/core/TableContainer';
import TableHead from '@material-ui/core/TableHead';
import TableRow from '@material-ui/core/TableRow';
import Paper from '@material-ui/core/Paper';
import { makeStyles } from '@material-ui/core';
import { sizing } from '@material-ui/system';


export class FormGeneDetails extends Component {
    continue = e => {
        e.preventDefault()
        this.props.nextStep()
    }

    sequenceChosen = e => {
        e.preventDefault()
        this.props.callAPI_spliceVariants()
        this.props.chooseGeneDetails()
    }


    render() {
        const { values, handleChange} = this.props
        // Function for dislaying after gene and species chosen 
        const renderVariants = () => {
            if(values.loading_spliceVariantsInfo) {
                return <h2>Loading splice variant information ...</h2>
            }
            if(values.loaded_spliceVariantsInfo) {

                const num_spliceVariants = values.spliceVariants_apiResponse.data.length 
                const array_spliceVariants = Array.from({length: num_spliceVariants}, (_,i) => i + 1)

                const options = []
                for (let i = 1; i <= num_spliceVariants; i += 1) {options.push(i);}

                return(
                    <React.Fragment>
                        <br/>    
                        <br/>    
                        <h1>Splice variants</h1>
                        <Grid container spacing={1} >
                            <Grid container item xs={12} direction='column'  alignItems = 'center' justify = 'center'>
                                <Box width='20' display="flex" flexDirection="column">
                                    <TableContainer component={Paper} align="center">
                                        <Table size="small" aria-label="a dense table">
                                            <TableHead>
                                                <TableRow>
                                                    <TableCell align="center">TranscriptNum</TableCell>
                                                    <TableCell align="center">TranscriptID</TableCell>
                                                    <TableCell align="center">Name</TableCell>
                                                    <TableCell align="center">Assembly</TableCell>
                                                    <TableCell align="center">Type</TableCell>
                                                    <TableCell align="center">AA_Length</TableCell>
                                                    <TableCell align="center">Is_Canonical</TableCell>
                                                </TableRow>
                                            </TableHead>
                                            {values.spliceVariants_apiResponse.data.map((row) => (
                                                <TableRow key={row.name}>
                                                    <TableCell align="center">{row[1]}</TableCell>
                                                    <TableCell align="center">{row.[2]}</TableCell>
                                                    <TableCell align="center">{row.[3]}</TableCell>
                                                    <TableCell align="center">{row.[4]}</TableCell>
                                                    <TableCell align="center">{row.[5]}</TableCell>
                                                    <TableCell align="center">{row.[6]}</TableCell>
                                                    <TableCell align="center">{row.[7]}</TableCell>
                                                </TableRow>
                                                ))}
                                            <TableBody>
                                            </TableBody>
                                        </Table>
                                    </TableContainer>
                                </Box>
                            </Grid>
                        </Grid>

                        <br/>
                        <br/>
                        <br/>    
                        <h1>Select variant and search sequence</h1>
                        <Box display="flex" flexDirection="row" p={1} m={1} style={{ width: '100%' }} alignItems="center" justifyContent="center">
                            <Box p={1}>
                                <InputLabel htmlFor="age-native-simple">Splice Variant</InputLabel>
                                <Select native onChange = {handleChange('spliceVariant')} inputinputProps = {{ style: {textAlign: 'center'}} }>
                                    {options.map(option => (
                                        <option key={option} value = {option}>
                                            {option} 
                                        </option>
                                    ))}
                                </Select>
                            </Box>
                            <Box p={1}>
                                <InputLabel htmlFor="age-native-simple">Search seq</InputLabel>
                                <Select native onChange = {handleChange('searchSeq')} inputinputProps = {{ style: {textAlign: 'center'}} }>
                                    <option value={'CDS'}>CDS</option>
                                    <option value={'cDNA'}>cDNA</option>
                                    ))}
                                </Select>
                            </Box>
                            <Box p={1}>
                                <Button
                                    label = "Continue"
                                    variant = "contained"
                                    style = {styles.button}
                                    color = "primary"
                                    onClick = {this.continue}
                                >
                                Continue</Button>
                            </Box>
                            <br/>
                        </Box>
                    </React.Fragment>
                )
            }
        }

        // this.props.values
        // Wrap everything in MuiThemeProvider
        // Change AppBar position later ... will not always display properly .. 
        return (
            <MuiThemeProvider>
                <React.Fragment>
                    <AppBar position = "sticky"> 
                        <Toolbar>
                            <Typography>
                                Automated sesRNA search
                            </Typography>
                        </Toolbar>
                    </AppBar>
                    <br/>
                    <h1>Enter Gene Details</h1>
                    <TextField 
                        placeholder = "Enter Species Name"
                        // floatingLabelText = "Species Name"
                        onChange = {handleChange('species')}
                        defaultValue={values.species}
                        style = {styles.textbox}
                    />
                    <br/>
                    <TextField 
                        placeholder = "Enter gene name"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('gene')}
                        defaultValue={values.gene}
                        style = {styles.textbox}
                    />
                    <br/>
                    <Button
                        label = "See availible splice variant"
                        variant = "contained"
                        style = {styles.button}
                        color = "secondary"
                        onClick = {this.sequenceChosen}
                    >
                    See availible splice variants
                    </Button>
                    {renderVariants()}
                </React.Fragment>
            </MuiThemeProvider>
        )
    }
}

const styles = {
    appbar: {
        margin: 20
    },
    textbox: {
        margin: 5
    },
    button: {
        margin: 15
    }, 
}


export default FormGeneDetails
