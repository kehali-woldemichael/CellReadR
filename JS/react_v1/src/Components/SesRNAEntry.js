import React, { Component } from 'react'
import FormGeneDetails from './FormGeneDetails'
import FormSesRNApar from './FormSesRNApar.js'
import EnterSequence from './EnterSequence'
import DisplaysesRNAs from './DisplaysesRNAs'

// For rendering json from node Express server 
import { JsonTable } from "react-json-to-table"
import {
    Card,
    CardBody,
    CardFooter,
    CardHeader,
    CardTitle,
    Col,
    Row,
  } from 'reactstrap';

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

        sesRNAs_apiResponse: "Default",
        loadingTable_sesRNAs: false, 
        loadedTable_sesRNAs: false
    }

    // Proceed to next step 
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
        this.setState({ [input]: e.target.value })
    }


    // When called ... fetches JSON of sesRNAs from node.js Express server 
    // Then changes loading state to false 
    callAPI_sesRNAs = async() => {
        this.setState({ loadingTable_sesRNAs: true } )

        const url = "http://localhost:9000"
        // const url = "https://api.randomuser.me/"

        const response = await fetch(url)
        const data = await response.json()

        this.setState({ sesRNAs_apiResponse: data } )

        this.setState({ loadingTable_sesRNAs: false } )
        this.setState({ loadedTable_sesRNAs: true } )
    }

    render_jsonTable = () => {
        if (this.state.sesRNAs_apiResponse && this.state.loadedTable_sesRNAs){
            var apiStatus = this.state.loadedTable_sesRNAs
            var apiData = this.state.sesRNAs_apiResponse
            if (apiStatus === true && apiData.length === 1){
                apiData = apiData[0]
                return(
                    <div>
                        <span className="text-success">Rows: <span className="text-primary"><b>{apiData.rows}</b></span></span>
                        <br/>
                        <span className="text-success">Columns: <span className="text-primary"><b>{apiData.cols}</b></span></span>
                        <br/>
                        <span className="text-success">Column Names: 
                            <span className="text-primary"><b>{this.renderColumnNames(apiData.columns)}</b></span>
                        </span>
                        <hr/>
                        <Card>
                            <CardHeader>All Rows</CardHeader>
                            <CardBody className="mb-1" style={{height:'400px', overflowY: "auto", overflow: "-moz-scrollbars-horizontal"}}>
                                <JsonTable rows={apiData.rowData} columns={apiData.columns} />
                            </CardBody>
                        </Card>
                    </div>
                )
            } else {
                return(
                    <div>
                        <span className="text-danger">{apiData[0].message}</span>
                    </div>
                )
            }
        }
    }

    render() {
        // Pull step out of state 
        const { step } = this.state
        // Pull out fields 
        const {species, gene, 
            spliceVariant, searchSeq,
            seqDirection, len_sesRNA, minTGG, maxStop, minGC, maxGC, 
            dist_cTGG, dist_stop_cTGG, choice_ATG, sesRNAs_apiResponse, loadingTable_sesRNAs, loadedTable_sesRNAs } = this.state 
        // To pass values to each component 
        const values = {species, gene, 
            spliceVariant, searchSeq,
            seqDirection, len_sesRNA, minTGG, maxStop, minGC, maxGC, 
            dist_cTGG, dist_stop_cTGG, choice_ATG, sesRNAs_apiResponse, loadingTable_sesRNAs, loadedTable_sesRNAs }

        switch(step) {
            case 1: 
                return (
                    <EnterSequence
                        // Props to access next step, handle change, values
                        nextStep = {this.nextStep}
                        handleChange = {this.handleChange}
                        values = {values}
                    />
                )
            // Enter gene information 
            case 2: 
                return (
                    <FormGeneDetails
                        // Props to access next step, handle change, values
                        nextStep = {this.nextStep}
                        prevStep = {this.prevStep}
                        handleChange = {this.handleChange}
                        values = {values}
                    />
                )
            // Enter sesRNA parameters
            case 3: 
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
            case 4: 
                return (
                    // Props to access confirm or previous step 
                    <DisplaysesRNAs
                        // nextStep = {this.nextStep}
                        prevStep = {this.prevStep}
                        handleChange = {this.handleChange}
                        values = {values}

                        // Calls API for returning sesRNAs 
                        callAPI_sesRNAs = {this.callAPI_sesRNAs}
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