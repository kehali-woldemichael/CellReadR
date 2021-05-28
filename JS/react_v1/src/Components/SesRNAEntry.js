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
        choice_ATG: '',         

        minTGG: '',         
        maxStop: '',         
        minGC: '',         
        maxGC: '',         
        dist_cTGG: '',         
        dist_stop_cTGG: '',         

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
            choice_ATG: 'None',         
            minTGG: '2',         
            maxStop: '0',         
            minGC: '30',         
            maxGC: '70',         
            dist_cTGG: '10',         
            dist_stop_cTGG: '20'
        })
    }
    autofill_medium = () => {
        this.setState({
            seqDirection: 'Both',         
            len_sesRNA: '204',         
            choice_ATG: 'All upstream',
            minTGG: '1',         
            maxStop: '1',         
            minGC: '30',         
            maxGC: '75',         
            dist_cTGG: '20',         
            dist_stop_cTGG: '10'         
        })
    }
    autofill_permissive = () => {
        this.setState({
            seqDirection: 'Both',         
            len_sesRNA: '204',         
            choice_ATG: 'Upstream central TGG',
            minTGG: '1',         
            maxStop: '2',         
            minGC: '30',         
            maxGC: '75',         
            dist_cTGG: '30',         
            dist_stop_cTGG: '5'         
        })
    }


    // When called ... fetches JSON of sesRNAs from node.js Express server 
    // Then changes loading state to false 
    callAPI_sesRNAs = async() => {
        // Changing state in order display loading screen 
        this.setState({ loading_sesRNAs: true } )
        // POST request using axios with async/await
        const sesRNA_parameters = {   "species" : this.state.species, 
                                "gene" : this.state.gene,
                                "spliceVariant" : this.state.spliceVariant,
                                "variantTable" : this.state.spliceVariants_apiResponse, 
                                "searchSeq" : this.state.searchSeq,
                                "seqDirection" : this.state.seqDirection,
                                "len_sesRNA" : this.state.len_sesRNA,
                                "minTGG" : this.state.minTGG,
                                "maxStop" : this.state.maxStop, 
                                "choice_ATG" : this.state.choice_ATG,
                                "minGC" : this.state.minGC, 
                                "maxGC" : this.state.maxGC, 
                                "dist_cTGG" : this.state.dist_cTGG,
                                "dist_stop_cTGG" : this.state.dist_stop_cTGG}
        console.log(sesRNA_parameters)

        const url = "http://localhost:9000/POST_sesRNAs"
        const response = await axios.post(url, sesRNA_parameters)
        console.log(response)

        this.setState({ sesRNAs_apiResponse: response } )
        this.setState({ loading_sesRNAs: false } )
        this.setState({ loaded_sesRNAs: true } )
    }

    callAPI_spliceVariants = async() => {
        // Changing state in order display loading screen 
        this.setState({ loading_spliceVariantsInfo: true } )
        // POST request using axios with async/await
        const geneDetails = { "species": this.state.species, "gene": this.state.gene}
        console.log(geneDetails)

        const url = "http://localhost:9000/POST_spliceVariants"
        const response = await axios.post(url, geneDetails)
        console.log(response)

        this.setState({ spliceVariants_apiResponse: response } )
        this.setState({ loading_spliceVariantsInfo: false } )
        this.setState({ loaded_spliceVariantsInfo: true } )
    }

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
            seqDirection, len_sesRNA, minTGG, maxStop, choiceATG, 
            minGC, maxGC, dist_cTGG, dist_stop_cTGG, 
            chosenGene, loading_spliceVariantsInfo, loaded_spliceVariantsInfo, spliceVariants_apiResponse,
            loading_sesRNAs, loaded_sesRNAs, sesRNAs_apiResponse} = this.state 
        // To pass values to each component 
        const values = {species, gene, 
            spliceVariant, searchSeq,
            seqDirection, len_sesRNA, minTGG, maxStop, choiceATG,
            minGC, maxGC, dist_cTGG, dist_stop_cTGG, 
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

                        callAPI_sesRNAs = {this.callAPI_sesRNAs}
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