import React, { Component } from 'react'
import FormGeneDetails from './FormGeneDetails'
import FormSesRNApar from './FormSesRNApar.js'
import ConfirmSesRNApar from './ConfirmSesRNApar.js'

export class SesRNAEntry extends Component {
    state = {
        step: 1,

        species: '',
        gene: '',

        spliceVariant: '',         
        searchSeq: '',         

        seqDirection: '',         
        len_sesRNA: '',         
        minTGG: '',         
        maxStop: '',         
        minGC: '',         
        maxGC: '',         
        dist_cTGG: '',         
        dist_stop_cTGG: '',         
        choice_ATG: '',         
    }

    // Proceed to next step 
    // Pulling out step from state and putting that into a variable 

    nextStep = () => {
        const { step } = this.state;
        this.setState({
        step: step + 1
        })
    }

    prevStep = () => {
        const { step } = this.state; 
        this.setState({
            step: step - 1
        }) 
    }

    //Handle field change 
    handleChange = input => e => {
        this.setState({ [input]: e.target.value });
    };

    render() {
        // Pull step out of state 
        const { step } = this.state
        // Pull out fields 
        const {species, gene, 
            spliceVariant, searchSeq,
            seqDirection, len_sesRNA, minTGG, maxStop, minGC, maxGC, 
            dist_cTGG, dist_stop_cTGG, choice_ATG } = this.state 
        // To pass values to each component 
        const values = {species, gene, 
            spliceVariant, searchSeq,
            seqDirection, len_sesRNA, minTGG, maxStop, minGC, maxGC, 
            dist_cTGG, dist_stop_cTGG, choice_ATG }

        switch(step) {
            // Enter gene information 
            case 1: 
                return (
                    <FormGeneDetails
                        // Props to access next step, handle change, values
                        nextStep = {this.nextStep}
                        handleChange = {this.handleChange}
                        values = {values}
                    />
                )
            // Enter sesRNA parameters
            case 2: 
                return (
                    <FormSesRNApar
                        // Props to access next step, previous step, handle change, values
                        nextStep = {this.nextStep}
                        prevStep = {this.prevStep}
                        handleChange = {this.handleChange}
                        values = {values}
                    />
                )
            // Confirm sesRNA parameters
            case 3: 
                return (
                    // Props to access confirm or previous step 
                    <ConfirmSesRNApar
                        // nextStep = {this.nextStep}
                        prevStep = {this.prevStep}
                        handleChange = {this.handleChange}
                        values = {values}
                    />
                )
        }

        return (
            <div>
                
            </div>
        )
    }
}

export default SesRNAEntry