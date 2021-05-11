import os
import pandas as pd
from kCellReadR.sequence import *

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
    def __init__(self, typeSeq, isoform, length,
                 num_inF_TGG, num_inF_Stop, num_inF_ATG,
                 minGC, maxGC, nearCenter, fromStop):
        self.typeSeq = typeSeq
        self.isoform = isoform
        self.length = length
        self.num_inF_TGG = num_inF_TGG
        self.num_inF_Stop = num_inF_Stop
        self.num_inF_ATG = num_inF_ATG
        self.minGC = minGC
        self.maxGC = maxGC
        self.nearCenter = nearCenter
        self.fromStop = fromStop

    def print_parameters(self):
        print('[instance attributes]')
        for attribute, value in self.__dict__.items():
            print(attribute, '=', value)

def return_inCDS(sesRNAs, CDS, isoForm, typeSes):
    """Return sesRNAs that are in CDS"""
    temp_cds_sesRNAs = []

    for sesRNA in sesRNAs:
        if typeSes == 'Reverse':
            if 0 != CDS[isoForm].seq.count(sesRNA.reverse_complement()):
                temp_cds_sesRNAs.append(sesRNA)
        elif typeSes == 'Complement':
            if 0 != CDS[isoForm].seq.count(sesRNA.complement()):
                temp_cds_sesRNAs.append(sesRNA)

    return temp_cds_sesRNAs

def generate_all_sesRNAs(rC_Seq, C_Seq, searchSequence, parameters):
    """Functoion for generating sesRNAs for both complement and reverse
    :rC_Seq: TODO
    :C_Seq: TODO
    :parameters: TODO
    :returns: TODO
    """

    if parameters.typeSeq == 'Reverse':
        rC_sesRNAs, rC_sequenceMetrics, rC_sesRNA_objs = \
            generate_sesRNAs_multiExon(rC_Seq, searchSequence, parameters)
        return rC_sesRNAs, rC_sequenceMetrics, rC_sesRNA_objs
    elif parameters.typeSeq == 'Complement':
        C_sesRNAs, C_sequenceMetrics, C_sesRNA_objs = \
            generate_sesRNAs_multiExon(C_Seq, searchSequence, parameters)
        return C_sesRNAs, C_sequenceMetrics, C_sesRNA_objs
    elif parameters.typeSeq == 'Both':
        all_sesRNAs = []
        all_sesRNA_objs = []

        parameters.typeSeq = 'Reverse'
        rC_sesRNAs, rC_sequenceMetrics, rC_sesRNA_objs = \
            generate_sesRNAs_multiExon(rC_Seq, searchSequence, parameters)
        all_sesRNAs.extend(rC_sesRNAs)
        all_sesRNA_objs.extend(rC_sesRNA_objs)

        parameters.typeSeq = 'Complement'
        C_sesRNAs, C_sequenceMetrics, C_sesRNA_objs = \
            generate_sesRNAs_multiExon(C_Seq, searchSequence, parameters)
        all_sesRNAs.extend(C_sesRNAs)
        all_sesRNA_objs.extend(C_sesRNA_objs)

        all_sequenceMetrics = \
            pd.concat([rC_sequenceMetrics, C_sequenceMetrics])
        return all_sesRNAs, all_sequenceMetrics, all_sesRNA_objs


def generate_sesRNAs_multiExon(exon_records, searchSequence, parameters):
    """Returs sesRNAs for all exons"""
    # Empty list for storing
    all_sesRNAs = []
    # Create empty dataframe for storting cell metrics
    all_sequenceMetrics = pd.DataFrame()
    all_sesRNA_objs = []

    # Iterating through exons and generating sesRNAs and sequence metrics
    for i in range(len(exon_records)):
        tempSeq = exon_records[i].seq
        temp_sesRNAs, temp_sequenceMetrics, temp_sesRNA_objs = \
            generate_sesRNA(tempSeq, searchSequence, parameters, (i+1))

        # Addiing sesRNAs and cell metrics for this exon
        all_sesRNAs.extend(temp_sesRNAs)
        all_sesRNA_objs.extend(temp_sesRNA_objs)
        all_sequenceMetrics = all_sequenceMetrics.append(temp_sequenceMetrics)

        # Printing number of passed sequences for current exon
        print(len(temp_sesRNAs))

    # Creating column for number of sequence ...
    # Then moving that column to the front
    all_sequenceMetrics['SeqNumber'] = list(range(1,(len(all_sequenceMetrics))+1))
    col = all_sequenceMetrics.pop("SeqNumber")
    all_sequenceMetrics.insert(0, col.name, col)

    # Return final output
    return all_sesRNAs, all_sequenceMetrics, all_sesRNA_objs

def generate_sesRNA(sequence, searchSequence, parameters, exonNumber):
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
    num_inF_ATGs = []
    num_inF_Stops = []

    sesRNA_objs = []

    while(start <= (len(sequence) - length)):
        # Defining current sub-sequence to process
        subsequence = sequence[start:(start+length)]

        # GC content
        gcContent = metric_gcContent(subsequence)*100
        num_inF_TGG, num_inF_ATG, num_inF_Stop, indices_inF_TGG, \
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
        if parameters.num_inF_ATG == 'None':
            cond3 = num_inF_ATG == 0
        elif parameters.num_inF_ATG == "All upstream":
            if num_inF_TGG != 0 and num_inF_ATG != 0:
                # Make sure all in frame ATG's upstream of all in frame TGG's
                cond3 = (min(indices_inF_TGG) > max(indices_inF_ATG))
            else:
                cond3 = num_inF_ATG == 0
        elif parameters.num_inF_ATG == "Upstream central":
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

        # Checking if any central array is more than 10 away
        # from an in frame stop
        if num_inF_Stop != 0 and num_inF_TGG != 0:
            if num_inF_TGG > 1:
                cond7 = any((min(abs(indices_inF_Stop - i)) >= parameters.fromStop)
                            for i in centralTGGs)
            else:
                cond7 = any((min(abs(indices_inF_Stop - i)) >= parameters.fromStop)
                            for i in [centralTGGs[0]])
        else:
            cond7 = True

        if(cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7):
            # Only include if in region of gene (currently in CDS)
            if check_inCDS(subsequence, searchSequence,
                           parameters.isoform, parameters.typeSeq):

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
                num_inF_ATGs.append(num_inF_ATG)
                num_inF_Stops.append(num_inF_Stop)

                sesRNA_objs.append(sesRNA(subsequence, start, start+length,
                                          indices_inF_TGG[0], central_inF_TGG,
                                          centralTGGs[1], num_inF_TGGs,
                                          num_inF_ATGs, num_inF_Stops,
                                          gcContent))


        total += 1
        # Updating start index
        start += 1



    allMetrics = {'TypeSeq': parameters.typeSeq, 'Exon':exonNumber, 'StartSeq':startSeq, 'StopSeq':stopSeq,
                  'firstTGG': first_TGGs, 'centralTGG': most_centralTGGs,
                  'second_cTGG': second_centralTGGs,
                  'numTGG':num_inF_TGGs, 'numATG':num_inF_ATGs,
                  'numStop':num_inF_Stops, 'gcContent':gcContents}
    sequenceMetrics = pd.DataFrame(allMetrics)

    return sesSeq, sequenceMetrics, sesRNA_objs
