import React, { Component } from 'react'
import {MuiThemeProvider} from "@material-ui/core/styles"

import { TextField } from '@material-ui/core'

import { Button } from '@material-ui/core'
import { ButtonGroup } from '@material-ui/core'

import AppBar from '@material-ui/core/AppBar'
import { Toolbar } from '@material-ui/core'
import { Typography } from '@material-ui/core'
 
import { Input } from '@material-ui/core';
// 
import { FormControl } from '@material-ui/core';
// import { MenuIcon } from '@material-ui/icons/Menu'
// import { IconButton } from '@material-ui/core'
import { InputLabel } from '@material-ui/core'

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
                                Automated sesRNA search
                            </Typography>
                        </Toolbar>
                    </AppBar>
                    <br/>
                    <h1>Choose sesRNA parameters</h1>
                    <h2>Autofill parameters</h2>
                    <ButtonGroup>
                        <Button
                            label = "Low"
                            variant = "contained"
                            style = {styles.button}
                            color = "primary"
                            // onClick = {this.back}
                        >Low</Button>
                        <Button
                            label = "Medium"
                            variant = "contained"
                            style = {styles.button}
                            color = "secondary"
                            // onClick = {this.back}
                        >Medium</Button>
                        <Button
                            label = "Strict"
                            variant = "contained"
                            style = {styles.button}
                            color = "secondary"
                            // onClick = {this.back}
                        >Strict</Button>
                    </ButtonGroup>
                    <br/>
                    <br/>
                    <h2>Manually enter sesRNA parameters</h2>
                    <h3>Sequence Search</h3>
                    <InputLabel>Splice Variant</InputLabel>
                    <TextField 
                        // InputLabelProps= {{ style: {textAlign: 'left'}}}
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "1, 2, ..."
                        // floatingLabelText = "Splice Variant"
                        onChange = {handleChange('spliceVariant')}
                        defaultValue={values.spliceVariant}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Search sequence</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "CDS, cDNA"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('searchSeq')}
                        defaultValue={values.searchSeq}
                        style = {styles.textbox}
                    />
                    <br/>
                    <br/>
                    <h3>sesRNA Feature Parameters</h3>
                    <br/>
                    <InputLabel>Sequence Direction</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "Both, Reverse Complement, Complement"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('seqDirection')}
                        defaultValue={values.seqDirection}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Length of sesRNA</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "200-300 bp"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('len_sesRNA')}
                        defaultValue={values.len_sesRNA}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Minimum number of TGGs</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "1, 2, ..."
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('minTGG')}
                        defaultValue={values.minTGG}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Maximum number of Stop Codons</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "0, 1, 2"
                        onChange = {handleChange('maxStop')}
                        defaultValue={values.maxStop}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Minimum GC content</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "30%+"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('minGC')}
                        defaultValue={values.minGC}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Maximum GC content</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "<75%"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('maxGC')}
                        defaultValue={values.maxGC}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Central TGG minimum distance from center</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "<20 bp"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('dist_cTGG')}
                        defaultValue={values.dist_cTGG}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>Minimum distance of central TGG from Stop codon</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = ">5 bp"
                        //InputLabelProps = "Gene Name"
                        onChange = {handleChange('dist_stop_cTGG')}
                        defaultValue={values.dist_stop_cTGG}
                        style = {styles.textbox}
                    />
                    <br/>
                    <InputLabel>In frame ATGs</InputLabel>
                    <TextField 
                        inputProps = {{ style: {textAlign: 'center'} }}
                        placeholder = "None, All upstream, One upstream"
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
                            label = "Continue"
                            variant = "contained"
                            style = {styles.button}
                            color = "secondary"
                            onClick = {this.continue}
                        >Continue</Button>
                    </ButtonGroup>
                </React.Fragment>
            </MuiThemeProvider>
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
        width: '10%',
    },
    button: {
        margin: 15
    }
}


export default FormSesRNApar