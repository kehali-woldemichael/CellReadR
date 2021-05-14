import React, { Component } from 'react'
import {MuiThemeProvider} from "@material-ui/core/styles"

import { TextField } from '@material-ui/core'

import { Button } from '@material-ui/core'

import AppBar from '@material-ui/core/AppBar'
import { Toolbar } from '@material-ui/core'
import { Typography } from '@material-ui/core'
 
// import { MenuIcon } from '@material-ui/icons/Menu'
// import { IconButton } from '@material-ui/core'
// import { InputLabel } from '@material-ui/core'

export class FormSesRNApar extends Component {
    continue = e => {
        e.preventDefault()
        this.props.nextStep()
    }
    back = e => {
        e.preventDefault()
        this.props.prevStep();
    }

    render() {
        const { values, handleChange} = this.props
        // this.props.values
        // Wrap everything in MuiThemeProvider
        // Change AppBar position later ... will not always display properly .. 
        return (
            <MuiThemeProvider>
                <React.Fragment>
                    <AppBar position = "sticky"> 
                        <Toolbar>
                            <Typography>
                                Enter sesRNA parameters
                            </Typography>
                        </Toolbar>
                    </AppBar>
                    <br/>
                    <TextField 
                        placeholder = "Enter splice variant"
                        // floatingLabelText = "Splice Variant"
                        onChange = {handleChange('spliceVariant')}
                        defaultValue={values.spliceVariant}
                        style = {styles.textbox}
                    />
                    <br/>
                    <TextField 
                        placeholder = "Enter Search Sequence (CDS, CDNA)"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('searchSeq')}
                        defaultValue={values.searchSeq}
                        style = {styles.textbox}
                    />
                    <br/>
                    <Button
                        label = "Back"
                        variant = "contained"
                        style = {styles.button}
                        color = "primary"
                        onClick = {this.back}
                    >Back</Button>
                    <Button
                        label = "Continue"
                        variant = "contained"
                        style = {styles.button}
                        color = "secondary"
                        onClick = {this.continue}
                    >Continue</Button>
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
    }
}

export default FormSesRNApar
