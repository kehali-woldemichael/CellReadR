import React, { Component } from 'react'

export class sesRNA_Entry extends Component {
    state = {
        step: 1,
        species: '',
        gene: '',
        
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
        this.setState({[input]: e.target.value})
    }

    render() {
        // Pull step out of state 
        const { step } = this.state
        const {species, gene} = this.state 
        const values = {species, gene} 

        switchState(step) {
            case 1: 
                return (
                    <FormUserDetails
                        nextStep = {this.nextStep}
                        handleChange = {this.handleChange}
                        values = {values}
                    />
                )
            case 2:
                return <h1></h1>
        }

        return (
            <div>
                
            </div>
        )
    }
}

export default sesRNA_Entry