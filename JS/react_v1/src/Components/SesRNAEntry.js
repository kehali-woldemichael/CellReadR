import React, { Component } from 'react'
import axios from "axios"; // for HTTP requests

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

        species: 'Rat',
        gene: 'Fezf2',

        spliceVariant: '1',         

        searchSeq: 'CDS',         

        seqDirection: '',         
        len_sesRNA: '',         
        minTGG: '',         
        maxStop: '',         
        minGC: '',         
        maxGC: '',         
        dist_cTGG: '',         
        dist_stop_cTGG: '',         
        choice_ATG: '',         

        chosenGene: false,
        loading_spliceVariantsInfo : false,
        loaded_spliceVariantsInfo : false,
        spliceVariants_apiResponse: "Default",


        loading_sesRNAs: false, 
        loaded_sesRNAs: false,
        sesRNAs_apiResponse: "Default"
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

    chooseGene = () => {
        this.setState({ chosenGene: true })
    }

    autofill_strict = () => {
        this.setState({
            seqDirection: 'Both',         
            len_sesRNA: '204',         
            minTGG: '2',         
            maxStop: '0',         
            minGC: '30',         
            maxGC: '70',         
            dist_cTGG: '10',         
            dist_stop_cTGG: '20',         
            choice_ATG: 'None'         
        })
    }
    autofill_medium = () => {
        this.setState({
            seqDirection: 'Both',         
            len_sesRNA: '204',         
            minTGG: '1',         
            maxStop: '1',         
            minGC: '30',         
            maxGC: '75',         
            dist_cTGG: '20',         
            dist_stop_cTGG: '10',         
            choice_ATG: 'All upstream'         
        })
    }
    autofill_permissive = () => {
        this.setState({
            seqDirection: 'Both',         
            len_sesRNA: '204',         
            minTGG: '1',         
            maxStop: '2',         
            minGC: '30',         
            maxGC: '75',         
            dist_cTGG: '30',         
            dist_stop_cTGG: '5',         
            choice_ATG: 'Upstream central TGG'         
        })
    }


    // When called ... fetches JSON of sesRNAs from node.js Express server 
    // Then changes loading state to false 
    callAPI_sesRNAs = async() => {
        this.setState({ loadingTable_sesRNAs: true } )

        const url = "http://localhost:9001"
        // const url = "https://api.randomuser.me/"

        const response = await fetch(url)
        const data = await response.json()

        this.setState({ sesRNAs_apiResponse: data } )

        this.setState({ loadingTable_sesRNAs: false } )
        this.setState({ loadedTable_sesRNAs: true } )
    }
    callAPI_spliceVariants = async() => {
        this.setState({ loading_spliceVariantsInfo: true } )

        const url = "http://localhost:9000/GET"

        const response = await fetch(url)
        const data = await response.json()

        this.setState({ spliceVariants_apiResponse: data } )

        this.setState({ loading_spliceVariantsInfo: false } )
        this.setState({ loaded_spliceVariantsInfo: true } )
    }

    send_geneDetails = async() => {
        const url = "http://localhost:9000/POST"
            // POST request using axios with async/await
            const geneDetails = { gene: this.state.gene };
            const response = await axios.post(url, geneDetails);
            // this.setState({ articleId: response.data.id });
    };


    render_sesRNAs = () => {
        if (this.state.sesRNAs_apiResponse && this.state.loadedTable_sesRNAs){
            var apiStatus = this.state.loaded_sesRNAs
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
            dist_cTGG, dist_stop_cTGG, choice_ATG, 
            chosenGene, loading_spliceVariantsInfo, loaded_spliceVariantsInfo, spliceVariants_apiResponse,
            loading_sesRNAs, loaded_sesRNAs, sesRNAs_apiResponse} = this.state 
        // To pass values to each component 
        const values = {species, gene, 
            spliceVariant, searchSeq,
            seqDirection, len_sesRNA, minTGG, maxStop, minGC, maxGC, 
            dist_cTGG, dist_stop_cTGG, choice_ATG, 
            chosenGene, loading_spliceVariantsInfo, loaded_spliceVariantsInfo, spliceVariants_apiResponse,
            loading_sesRNAs, loaded_sesRNAs, sesRNAs_apiResponse} 

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
                        values = {values}
                        nextStep = {this.nextStep}
                        prevStep = {this.prevStep}
                        chooseGeneDetails = {this.chooseGene}
                        handleChange = {this.handleChange}

                        callAPI_spliceVariants = {this.callAPI_spliceVariants}
                        send_geneDetails = {this.send_geneDetails}
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

                        autofill_strict = {this.autofill_strict}
                        autofill_medium = {this.autofill_medium}
                        autofill_permissive = {this.autofill_permissive}
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