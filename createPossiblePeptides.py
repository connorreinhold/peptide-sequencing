##WRITTEN BY CONNOR REINHOLD
import random
import numpy
import time

start_time=time.time()
def calculateWeight(line, containsWater, numProtons):
    proton_mass = 1.00728
    water_mass = 18.010565
    pm = {'G': 57.021464, 'D': 115.02694, 'A': 71.037114, 'Q': 128.05858, 'S': 87.032029,
          'K': 128.09496, 'P': 97.052764, 'E': 129.04259, 'V': 99.068414, 'M': 131.04048,
          'T': 101.04768, 'H': 137.05891, 'C': 103.00919, 'F': 147.06841, 'L': 113.08406,
          'R': 156.10111, 'I': 113.08406, 'N': 114.04293, 'Y': 163.06333, 'W': 186.07931,
          'm': 147.03539, 'U': 167.956}
    weight = 0
    for letter in line:
        weight += pm[letter]
    if containsWater:
        weight += water_mass
    weight += numProtons * proton_mass
    return weight


def createSubsetPeptideFile(peptide_file_name, fasta_file_name):
    with open(fasta_file_name) as fasta_file:
        sample_peptide_rate=0.000001  #draw 1 peptide for each whatever amino acids in the sequence
        peptide_file=open(peptide_file_name,"w")
        big_sequence=""

        def calculateSubset(sequence):
            if len(sequence) > 0:
                sequence = sequence.replace('\n', '')
                for i in range(0,len(sequence)-8):
                    random_number=random.random()
                    if random_number<sample_peptide_rate:
                        peptide_file.write(sequence[i:i + 9] + '\n')
                for i in range(0,len(sequence)-9):
                    random_number=random.random()
                    if random_number<sample_peptide_rate:
                        peptide_file.write(sequence[i:i + 10]+'\n')


        for line in fasta_file:
            if line[0]==">":
                calculateSubset(big_sequence)
                big_sequence=""

            else:
                big_sequence+=line

        calculateSubset(big_sequence)
        peptide_file.write(big_sequence[0:10])
        peptide_file.close()
    fasta_file.close()


def createReversedFile(reversed_fasta_file_name, fasta_file_name):
    with open(fasta_file_name) as fasta_file:
        reversed_fasta_file = open(reversed_fasta_file_name, "w")
        sequence=""
        def makeDecoy(decoy_sequence):
            if len(decoy_sequence) > 0:
                decoy_sequence = decoy_sequence.replace('\n', '')
                decoy_sequence = decoy_sequence[::-1]
                return '\n' + decoy_sequence
            else:
                return ""

        for line in fasta_file:
            if line[0]==">":
                if len(sequence)>0:
                    reversed_fasta_file.write(makeDecoy(sequence)+"\n")
                reversed_fasta_file.write(">reversed"+line[1:-1])
                sequence=""
            else:
                sequence+=line
        reversed_fasta_file.write(makeDecoy(sequence))
    reversed_fasta_file.close()
    fasta_file.close()

#nonmers and tenmers
def createCompletePeptideArrayFile(complete_peptide_list_name, fasta_file_name):
    with open(fasta_file_name) as fasta_file:
        complete_peptide_file=open(complete_peptide_list_name, "w")
        listPeptides=[]
        listWeights=[]
        sequence=""
        line_number=0
        for line in fasta_file:
            line_number+=1
            if line_number%2000==0:
                print(line_number, time.time()-start_time)
            if line[0]==">":
                if len(sequence)>0:
                    sequence=sequence.replace('\n','')
                    for nonindex in range(0,len(sequence)-8):
                        nonmer=sequence[nonindex:nonindex+9]
                        try:
                            listWeights.append(round(calculateWeight(nonmer, True, 2)/2,3))
                            listPeptides.append(nonmer)
                        except KeyError:
                            print("Could not recognize a character in " + nonmer)

                    for decaindex in range(0, len(sequence) - 9):
                        decamer = sequence[decaindex:decaindex + 10]
                        try:
                            listWeights.append(round(calculateWeight(decamer, True, 2) / 2, 3))
                            listPeptides.append(decamer)
                        except KeyError:
                            print("Could not recognize a character in "+decamer)
                    sequence=""
            else:
                sequence+=line
        sequence = sequence.replace('\n', '')
        for nonindex in range(0, len(sequence) - 8):
            nonmer = sequence[nonindex:nonindex + 9]
            try:
                listWeights.append(round(calculateWeight(nonmer, True, 2) / 2, 3))
                listPeptides.append(nonmer)
            except KeyError:
                print("Could not recognize a character in " + nonmer)
        for decaindex in range(0, len(sequence) - 9):
            decamer = sequence[decaindex:decaindex + 10]
            try:
                listWeights.append(round(calculateWeight(decamer, True, 2) / 2, 3))
                listPeptides.append(decamer)
            except KeyError:
                print("Could not recognize a character in " + decamer)
        listPeptides=numpy.array(listPeptides)
        listWeights=numpy.array(listWeights)
        print("Removing duplicates")
        listPeptides,dups=numpy.unique(listPeptides, return_index=True)
        listWeights=listWeights[dups]
        print("Finding indices...")
        inds=listWeights.argsort()
        print("Sorted weights - remembering indices...")
        listWeights=listWeights[inds]
        listPeptides=listPeptides[inds]
        print(len(listWeights))
        #listWeights=listWeights[0:len(listWeights)/2]
        #listPeptides = listPeptides[0:len(listPeptides) / 2]

        for index in range(0, len(listWeights)):
            if index==(len(listWeights)-1):
                complete_peptide_file.write(listPeptides[index] + "," + str(listWeights[index]))
            else:
                complete_peptide_file.write(listPeptides[index] + "," + str(listWeights[index]) + "\n")
    complete_peptide_file.close()
    fasta_file.close()


def main():
    start_time=time.time()
    fasta_file="human_proteome.fasta"

    #createCompletePeptideArrayFile("complete_peptide_file.txt", fasta_file)

    #createReversedFile("reversed_fasta_file.fasta", fasta_file)

    #createCompletePeptideArrayFile("complete_peptide_reversed_file.txt", "reversed_fasta_file.fasta")

    createSubsetPeptideFile("peptide_file.txt", fasta_file)

    print(time.time()-start_time)


if __name__ == '__main__':
    main()
