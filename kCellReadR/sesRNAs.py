import os
import pandas as pd
from kCellReadR.sequence import *
from kCellReadR.paths import *

# Class to store sesRNA object
class sesRNA():
    """Docustring for sesRNA class"""
    def __init__(self, seq, startSeq, stopSeq,
                 firstTGG, centralTGG, second_cTGG,
                 numTGG, numATG, numStop, gcContent):

        self.seq = seq
        self.startSeq = startSeq
        self.stopSeq = stopSeq
        self.firstTGG = firstTGG
        self.centralTGG = centralTGG
        self.second_cTGG = second_cTGG
        self.numTGG = numTGG
        self.numATG = numATG
        self.numStop = numStop
        self.gcContent = gcContent

class parameters_sesRNA:
    """Class for storing paramaters"""
    def __init__(self, species, gene, spliceVariant, length,
                 num_inF_TGG, num_inF_Stop, inF_ATG,
                 minGC, maxGC, nearCenter, fromStop):
        self.species = species
        self.gene = gene
        self.spliceVariant = spliceVariant

        self.length = length
        # Limit on number of in frame TGG's 
        self.num_inF_TGG = num_inF_TGG
        # Limit on number of in frame STOP codons (TAG, TGA, TAA) 
        self.num_inF_Stop = num_inF_Stop
        # Limit of number of in frame start codons (ATG) 
        self.inF_ATG = inF_ATG
        self.minGC = minGC
        self.maxGC = maxGC
        # Distance of most central TGG from center 
        self.nearCenter = nearCenter
        # Distance of most central TGG from any STOP
        self.fromStop = fromStop

    def print_parameters(self):
        print('[instance attributes]')
        for attribute, value in self.__dict__.items():
            print(attribute, '=', value)

def return_inSearchSeq(sesRNAs, searchSeq):
    """Return sesRNAs that are in given target sequence"""
    temp_cds_sesRNAs = []

    for sesRNA in sesRNAs:
        if 0 != searchSeq[0].seq.count(Seq(sesRNA).reverse_complement()):
            temp_cds_sesRNAs.append(sesRNA)

    return temp_cds_sesRNAs

def generate_all_sesRNAs(target, searchSequence, parameters, variantTable):
    """Function for generating sesRNAs for both complement and reverse
    :rC_Seq: reverse compelement of 
    :parameters: TODO
    :returns: TODO
    """
    rC_sesRNAs, rC_sequenceMetrics, rC_sesRNA_objs = \
        generate_sesRNAs_multiExon(target, searchSequence, parameters, variantTable)
    return rC_sesRNAs, rC_sequenceMetrics, rC_sesRNA_objs


def generate_sesRNAs_multiExon(exon_records, searchSequence, parameters, variantTable):
    """Returs sesRNAs for all exons"""
    # Empty list for storing
    all_sesRNAs = []
    # Create empty dataframe for storting cell metrics
    all_sequenceMetrics = pd.DataFrame()
    all_sesRNA_objs = []

    # Iterating through exons and generating sesRNAs and sequence metrics
    for i in range(len(exon_records)):
        tempSeq = exon_records[i].seq
        # Only generating sesRNAs if length of exon greater than length of sesRNA 
        if len(tempSeq) >= parameters.length:
            temp_sesRNAs, temp_sequenceMetrics, temp_sesRNA_objs = \
                generate_sesRNA(tempSeq, searchSequence, parameters, (i+1), variantTable)

            # Addiing sesRNAs and cell metrics for this exon
            all_sesRNAs.extend(temp_sesRNAs)
            all_sesRNA_objs.extend(temp_sesRNA_objs)
            all_sequenceMetrics = all_sequenceMetrics.append(temp_sequenceMetrics)

            i += 1
            print(i)
        # Skip if current exon too short to generate sesRNAs 
        else:
            i += 1

        # Printing number of passed sequences for current exon
        # print(len(temp_sesRNAs))

    # Creating column for number of sequence ...
    # Then moving that column to the front
    all_sequenceMetrics['sesNum'] = list(range(1,(len(all_sequenceMetrics))+1))
    col = all_sequenceMetrics.pop("sesNum")
    all_sequenceMetrics.insert(0, col.name, col)

    # Return final output
    return all_sesRNAs, all_sequenceMetrics, all_sesRNA_objs

def generate_sesRNA(sequence, searchSequence, parameters, exonNumber, variantTable):
    """Generates sesRNAs given sequence"""
    total = 0
    start = 0
    length = parameters.length
    center = length/2

    numTGG = []

    sesSeq = []
    startSeq = []
    stopSeq = []
    gcContents = []

    first_TGGs = []
    most_centralTGGs = []
    second_centralTGGs = []

    # For storing number of in frame TGG, ATG, and Stop codons
    num_inF_TGGs = []
    num_inF_TTGGs = []
    num_inF_TGGAs = []
    num_inF_TTGGAs = []

    num_inF_ATGs = []
    num_inF_Stops = []

    # If in exon object 
    exonTotal_Ratio = []
    exonProtein_Ratio = []
    cdsTotal_Ratio = []
    cdnaTotal_Ratio = []

    sesRNA_objs = []

    while(start <= (len(sequence) - length)):
        # Defining current sub-sequence to process
        subsequence = sequence[start:(start+length)]

        # GC content
        gcContent = metric_gcContent(subsequence)*100
        num_inF_TGG, num_inF_TTGG, num_inF_TGGA, num_inF_TTGGA, num_inF_ATG, num_inF_Stop, indices_inF_TGG, \
            indices_inF_ATG, indices_inF_Stop = \
            return_inFrame(subsequence, 'all')

        # Identifying most and second most central in frame TGG's
        if num_inF_TGG != 0:
            sorted_TGGs = list(
                (np.array(sorted(indices_inF_TGG - length/2, key = abs))
                 + (length/2)))
            central_inF_TGG = int(sorted_TGGs[0])
            if num_inF_TGG > 1:
                secondCentral_inF_TGG = int(sorted_TGGs[1])
                centralTGGs = [central_inF_TGG, secondCentral_inF_TGG]
            else:
                centralTGGs = [central_inF_TGG, 'NA']
                secondCentral_inF_TGG = 'NA'

        # Only proceed if passed
        cond1 = num_inF_Stop <= parameters.num_inF_Stop
        cond2 = num_inF_TGG >= parameters.num_inF_TGG

        # Checking ATGs ... either no in frame or only if upstream of all TGGs
        if parameters.inF_ATG == 'None':
            cond3 = num_inF_ATG == 0
        elif parameters.inF_ATG == "All upstream":
            if num_inF_TGG != 0 and num_inF_ATG != 0:
                # Make sure all in frame ATG's upstream of all in frame TGG's
                cond3 = (min(indices_inF_TGG) > max(indices_inF_ATG))
            else:
                cond3 = num_inF_ATG == 0
        elif parameters.inF_ATG == "Upstream central":
            if num_inF_TGG != 0 and num_inF_ATG != 0:
                # Make sure all in frame ATG's upstream of all in frame TGG's
                cond3 = (central_inF_TGG > max(indices_inF_ATG))
            else:
                cond3 = num_inF_ATG == 0

        # Condition to chekc gcContent
        cond4 = gcContent >= parameters.minGC
        cond5 = gcContent <= parameters.maxGC

        # Checking if TGG near center of subsequence
        cond6 = any(abs(x - center) < parameters.nearCenter for x in indices_inF_TGG)

        # Checking if central TGG is more than some bp away
        # from an in frame stop
        if num_inF_Stop != 0 and num_inF_TGG != 0:
            #if num_inF_TGG > 1:
                #cond7 = any((min(abs(indices_inF_Stop - i)) >= parameters.fromStop)
                            #for i in centralTGGs)
            #else:
                #cond7 = any((min(abs(indices_inF_Stop - i)) >= parameters.fromStop)
                            #for i in [centralTGGs[0]])
            cond7 = any((min(abs(indices_inF_Stop - i)) >= parameters.fromStop)
                            for i in [centralTGGs[0]])
        else:
            cond7 = True

        if(cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7):
            # Only include if in region of gene (currently in CDS)
            if check_inSearchSeq(subsequence, searchSequence):

                # Appending passesed subsequences
                sesSeq.append(subsequence)
                # Appending indices of start for ses
                # (relative to start of exon)
                startSeq.append(start)
                stopSeq.append(start+length)
                gcContents.append(gcContent)

                first_TGGs.append(indices_inF_TGG[0])
                most_centralTGGs.append(central_inF_TGG)
                second_centralTGGs.append(centralTGGs[1])

                # Appending number of in frame TGGs, ATGs, and Stop codons
                num_inF_TGGs.append(num_inF_TGG)
                num_inF_TTGGs.append(num_inF_TTGG)
                num_inF_TGGAs.append(num_inF_TGGA)
                num_inF_TTGGAs.append(num_inF_TTGGA)

                num_inF_ATGs.append(num_inF_ATG)
                num_inF_Stops.append(num_inF_Stop)

                tempTotal, tempProtein, tempCDS, tempCDNA = check_inExonVariants(subsequence, parameters.species, parameters.gene, variantTable)
                exonTotal_Ratio.append(tempTotal)
                exonProtein_Ratio.append(tempProtein)
                cdsTotal_Ratio.append(tempCDS)
                cdnaTotal_Ratio.append(tempCDNA)

                sesRNA_objs.append(sesRNA(subsequence, start, start+length,
                                          indices_inF_TGG[0], central_inF_TGG,
                                          centralTGGs[1], num_inF_TGGs, 
                                          num_inF_ATGs, num_inF_Stops,
                                          gcContent))


        total += 1
        # Updating start index
        start += 1



    allMetrics = {'SeqNum':exonNumber, 
                  "ExonFrac":exonTotal_Ratio, "ExonProtFrac": exonProtein_Ratio, 
                  "CDSFrac":cdsTotal_Ratio, "CDNAFrac":cdnaTotal_Ratio, 
                  'StartSeq':startSeq, 'StopSeq':stopSeq,
                  'firstTGG': first_TGGs, 'centralTGG': most_centralTGGs,
                  'second_cTGG': second_centralTGGs,
                  'numTGG':num_inF_TGGs, 
                  "numTTGG": num_inF_TTGGs, "numTGGA": num_inF_TGGAs, "numTTGGA": num_inF_TTGGAs, 
                  'numATG':num_inF_ATGs,
                  'numStop':num_inF_Stops, 'gcCont':gcContents}
    sequenceMetrics = pd.DataFrame(allMetrics)

    return sesSeq, sequenceMetrics, sesRNA_objs
