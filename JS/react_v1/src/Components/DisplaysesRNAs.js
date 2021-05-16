import React, { Component } from 'react'
import Grid from '@material-ui/core/Grid'
import {MuiThemeProvider} from "@material-ui/core/styles"

import { TextField } from '@material-ui/core'

import { Button } from '@material-ui/core'
import { ButtonGroup } from '@material-ui/core'
import Box from '@material-ui/core/Box'

import AppBar from '@material-ui/core/AppBar'
import { Toolbar } from '@material-ui/core'
import { Typography } from '@material-ui/core'

export class Display extends Component {
    back = e => {
        e.preventDefault()
        this.props.prevStep();
    }

    render() {
        const { values, handleChange} = this.props

        return (
            <Grid container spacing={2} >
                <AppBar position = "sticky"> 
                    <Toolbar>
                        <Typography>
                            Automated sesRNA search
                        </Typography>
                    </Toolbar>
                </AppBar>
                <br/>
                <Grid container item xs={2} direction='column'  alignItems = 'center' justify = 'center'>
                    <MuiThemeProvider>
                        <React.Fragment>
                            <h1>Sequence Search</h1>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Enter splice variant"
                                // floatingLabelText = "Splice Variant"
                                onChange = {handleChange('spliceVariant')}
                                defaultValue={values.spliceVariant}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Enter Search Sequence (CDS, CDNA)"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('searchSeq')}
                                defaultValue={values.searchSeq}
                                style = {styles.textbox}
                            />
                            <br/>
                            <h1>sesRNA Feature Parameters</h1>
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Sequence Direction"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('seqDirection')}
                                defaultValue={values.seqDirection}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Length of sesRNA"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('len_sesRNA')}
                                defaultValue={values.len_sesRNA}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Minimum number of TGGs"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('minTGG')}
                                defaultValue={values.minTGG}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Maximum number of Stop codons"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('maxStop')}
                                defaultValue={values.maxStop}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Minimum GC content"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('minGC')}
                                defaultValue={values.minGC}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Maximum GC content"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('maxGC')}
                                defaultValue={values.maxGC}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Central TGG minimum distance from center"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('dist_cTGG')}
                                defaultValue={values.dist_cTGG}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "Central TGG minimum distance from Stop codon"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('dist_stop_cTGG')}
                                defaultValue={values.dist_stop_cTGG}
                                style = {styles.textbox}
                            />
                            <br/>
                            <TextField 
                                inputProps = {{ style: {textAlign: 'center'} }}
                                placeholder = "How many in frame ATGs to allow"
                                //InputLabelProps = "Gene Name"
                                onChange = {handleChange('choice_ATG')}
                                defaultValue={values.choice_ATG}
                                style = {styles.textbox}
                            />
                            <br/>
                            <ButtonGroup>
                                <Button
                                    label = "Back"
                                    variant = "contained"
                                    style = {styles.button}
                                    color = "primary"
                                    onClick = {this.back}
                                >Back</Button>
                                <Button
                                    label = "Reload"
                                    variant = "contained"
                                    style = {styles.button}
                                    color = "secondary"
                                    // onClick = {this.back}
                                >Reload</Button>
                            </ButtonGroup>
                        </React.Fragment>
                    </MuiThemeProvider>
                </Grid>
                <Grid container item xs={6} direction="column" >
                    <MuiThemeProvider>
                        <h1>Table Candidate sesRNAs</h1>
                        <React.Fragment>
                            <ButtonGroup>
                                <Button
                                    label = "Compute higher order features"
                                    variant = "contained"
                                    style = {styles.button}
                                    color = "secondary"
                                    // onClick = {this.back}
                                >Toggle View</Button>
                                <Button
                                    label = "Compute higher order features"
                                    variant = "contained"
                                    style = {styles.button}
                                    color = "secondary"
                                    // onClick = {this.back}
                                >Compute higher order features</Button>
                                <Button
                                    label = "Download sesRNAs"
                                    variant = "contained"
                                    style = {styles.button}
                                    color = "secondary"
                                    // onClick = {this.back}
                                >Download sesRNAs</Button>
                            </ButtonGroup>
                        </React.Fragment>
                    </MuiThemeProvider>
                </Grid>
            </Grid>
        )
    }
}

// Adjusting the margin of main page components 
const styles = {
    appbar: {
        margin: 20
    },
    textbox: {
        margin: 5,
        width: '40%',
    },
    button: {
        margin: 10
    }
}

export default Display
