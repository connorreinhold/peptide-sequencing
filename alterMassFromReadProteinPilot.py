
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
        try:
            weight += pm[letter]
        except:
            print(letter)
    if containsWater:
        weight += water_mass
    weight += numProtons * proton_mass
    return weight

def main():
    fname = "ref_msf.txt"  # mascot formatted file
    fname1 = "protein_pilot_list.txt"  # peptide list

    filep1 = open(fname1, "r")



    f_out = open('ref.txt', 'w')

    firstLine = True
    for peptideLine in filep1:
        pepArgs = peptideLine.split()
        peptide = pepArgs[0]
        mz_tmp = pepArgs[1]
        scan = round(float(pepArgs[2]), 3)
        charge = int(pepArgs[3])
        amp_s = pepArgs[4]
        ret_time = pepArgs[5]
        # print(scan)

        foundPeptide=False
        peptideMass=calculateWeight(peptide,True,2)
        if charge==2:
            filep = open(fname, "r")
            for line_from in filep:
                if line_from[0:8] == "PEPMASS=":
                    mz_amp_msf = line_from.split()
                    scan_read = round(float(mz_amp_msf[3]), 3)
                    if scan == scan_read:
                        foundPeptide=True
                        if firstLine:
                            f_out.write("BEGIN IONS")
                            firstLine=False
                        else:
                            f_out.write("\n\nBEGIN IONS")
                        f_out.write("\nTITLE= "+peptide+ " ret_time "+ret_time)
                        f_out.write("\nPEPMASS= "+mz_tmp + " "+amp_s+" "+str(scan))
                        print(peptide)
                if foundPeptide:
                    if line_from[0:8] == "END IONS":
                        foundPeptide = False
                        f_out.write("\nEND IONS")
                    elif line_from[0:8] != "PEPMASS=":
                        fragment_mass=line_from.strip()
                        fragment_mass=float(fragment_mass[0:fragment_mass.find(' ')])
                        if abs(peptideMass/2-fragment_mass)>=2.5:
                            f_out.write('\n'+str(fragment_mass))
            filep.close()

    filep.close()
    filep1.close()
    f_out.close()


if __name__ == '__main__':
    main()
