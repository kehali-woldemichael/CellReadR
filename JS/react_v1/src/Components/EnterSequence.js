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

export class EnterSequence extends Component {
    continue = e => {
        e.preventDefault()
        this.props.nextStep()
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
                                Automated sesRNA search
                            </Typography>
                        </Toolbar>
                    </AppBar>
                    <h1>Choose method of sequence selection</h1>
                    <Button
                        label = "Continue"
                        variant = "contained"
                        style = {styles.button}
                        color = "secondary"
                        onClick = {this.continue}
                    >Use reference sequences</Button>
                    <Button
                        label = "Continue"
                        variant = "contained"
                        style = {styles.button}
                        color = "primary"
                        // onClick = {this.continue}
                    >Upload sequence</Button>
                    <Button
                        label = "Continue"
                        variant = "contained"
                        style = {styles.button}
                        color = "primary"
                        // onClick = {this.continue}
                    >Load previous trial</Button>
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

export default EnterSequence
