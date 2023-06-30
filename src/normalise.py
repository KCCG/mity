"""Normalise mity VCF."""

import sys
import gzip
import re
import logging
from .util import write_vcf
from scipy.stats import binom
from math import isinf
import numpy as np

def unchanged(List):
    # check that all numbers in the list are the same
    if len(set(List)) != 1:
        sys.exit("Multiple values in " + str(List))
    return ((List[0]))


def Sum(List):
    # convert to integer
    List = [float(i) for i in List]
    return (sum(List))


def weighted_mean(mean_list, count_list):
    mean_list = [float(i) for i in mean_list]
    count_list = [float(i) for i in count_list]

    numerator = [a * b for a, b in zip(mean_list, count_list)]
    numerator = sum(numerator)
    denominator = sum(count_list)
    return (numerator / denominator)


# this will split multiallelic lines in a vcf so that each allele is on a separate line
# variants is a list, each element is a line from the vcf
# eg variants = ['MT\t3103\t.\tCTACN\tCTAC,CTACT\t68191.9\t.\tAB=0,0;ABP=0,0;AC=2,4;AF=0.333333,0.666667;AN=6;;...,'MT\t3108\t.\tTT\tTC\t0\t.\tAB=0;ABP=0;AC=0;AF=0;AN=6;...]
def split_multi_allelic(variants):
    split_vcf = []

    # loop over each line in variants
    for line in variants:

        lineparts = line.split('\t')
        # print(lineparts)
        # lineparts = ['MT', '3103', '.', 'CTACN', 'CTAC,CTACT', '68191.9', '.', 'AB=0,0;ABP=0,0;AC=2,4;AF=0.333333,0.666667;AN=6;...]

        # The position of these fields should never change in a vcf:
        chromosome = lineparts[0]
        pos = lineparts[1]
        ID = lineparts[2]
        ref = lineparts[3]
        QUAL = lineparts[5]
        Filter = lineparts[6]

        alt = lineparts[4]
        # print(alt):
        # alt = CTAC,CTACT
        alt = alt.split(",")
        # print(alt):
        # alt = ['CTAC', 'CTACT']
        n_alleles = len(alt)

        ######################################
        ###### INFO field

        INFO = lineparts[7]
        INFO = INFO.split(";")
        # print(INFO):
        # INFO = ['AB=0,0', 'ABP=0,0', 'AC=2,4', 'AF=0.333333,0.666667', 'AN=6', ...]
        # get the abbreviated names in the info field
        INFO_names = [x.split("=")[0] for x in INFO]
        # print(INFO_names):
        # INFO_names = ['AB', 'ABP', 'AC', 'AF', 'AN'...]

        # INFO_to_keep are the INFO fields that you want to remain in the VCF.
        # Note that AF is replaced by VAF later, so don't need to keep this
        INFO_to_keep = ["DP", "QR", "QA", "SRF", "SRR", "SAF", "SAR", "MQM", "MQMR"]

        # get the index of where each of the INFO_to_keep are in the INFO
        # better way to do this?
        idx = []
        for keep_name in INFO_to_keep:
            counter = 0
            for names in INFO_names:
                # print(keep_name)
                if names == keep_name:
                    # print(names)
                    idx.append(counter)
                counter += 1
        idx = sorted(idx)
        # print(idx):
        # idx = [7, 15, 16, 26, 27, 34, 36, 37, 39]

        # pick out the info names, this should be the same as INFO_to_keep
        INFO_names = [INFO_names[i] for i in idx]
        # add the equal sign back in
        INFO_names = [x + "=" for x in INFO_names]
        # print(INFO_names):
        # INFO_names = ['DP=', 'MQM=', 'MQMR=', 'QA=', 'QR=', 'SAF=', 'SAR=', 'SRF=', 'SRR=']

        # get the values of the INFO
        INFO_values = [x.split("=")[1] for x in INFO]
        INFO_values = [INFO_values[i] for i in idx]
        # print(INFO_values)
        # INFO_values = ['7575', '70,60', '0', '75397,5623', '0', '1289,100', '1351,103', '0', '0']

        # split the values for each allele
        INFO_values = [x.split(",") for x in INFO_values]
        # print(INFO_values)
        # INFO_values = [['7575'], ['70', '60'], ['0'], ['75397', '5623'], ['0'], ['1289', '100'], ['1351', '103'], ['0'], ['0']]

        # get the number of values for each allele
        len_INFO_values = [len(x) for x in INFO_values]
        # print(len_INFO_values)
        # len_INFO_values = [1, 2, 1, 2, 1, 2, 2, 1, 1]

        # In INFO_values, there should only be two lengths of the sublists,
        # either length 1, if its a reference value, or length of the number of alleles
        if set(len_INFO_values) != set([1, n_alleles]):
            sys.exit(
                "the values in the info field must either be for the reference (length = one) or for the alternates (length = number of alleles)")

        # if there is only one value, i.e. just for the reference, repeat it the number of alleles times
        # so that it can be picked out later when each allele is put on a new line
        INFO_values = [x * n_alleles if len(x) == 1 else x for x in INFO_values]
        # info_values = [x*n_alleles if len(x)<n_alleles else x for x in info_values]
        # print(INFO_values)
        # INFO_values = [['7575', '7575'], ['70', '60'], ['0', '0'], ['75397', '5623'], ['0', '0'], ['1289', '100'], ['1351', '103'], ['0', '0'], ['0', '0']]

        ######################################
        ###### GENO/FORMAT field
        FORMAT_name = lineparts[8]
        FORMAT_name = FORMAT_name.split(":")
        logging.debug(FORMAT_name)
        # FORMAT_name = ['GT', 'DP', 'AD', 'RO', 'QR', 'AO', 'QA', 'GL']

        FORMAT_values = lineparts[9:]
        # print(FORMAT_values)
        #  = ['2/2:2915:0,0,88:0:0:0,88:0,2441:-251.884,-310.442,-249.634,-90.9111,-88.9622,0', '2/2:1936:0,0,67:0:0:0,67:0,1833:-174.582,-231.477,-174.582,-78.2678,-78.2678,0', '1/1:2724:0,2640,48:0:0:2640,48:75397,1349:-6804.65,-872.084,0,-6770.44,-762.912,-6680.73']

        FORMAT_values = [x.split(":") for x in FORMAT_values]
        # print(FORMAT_values)

        # remove GL as we dont really know what to do with it yet
        try:
            idx = FORMAT_name.index("GL")
            del FORMAT_name[idx]
            for x in FORMAT_values:
                del x[idx]
            # print(FORMAT_values)
            # FORMAT_values = [['2/2', '2915', '0,0,88', '0', '0', '0,88', '0,2441'], ['2/2', '1936', '0,0,67', '0', '0', '0,67', '0,1833'], ['1/1', '2724', '0,2640,48', '0', '0', '2640,48', '75397,1349']]
        except ValueError:
            pass
        n_samples = len(FORMAT_values)
        # print(n_samples)

        ##################################
        ####### Create new line

        for allele in range(0, n_alleles):

            # build up the new line
            new_line = [chromosome, pos, ID, ref, alt[allele], QUAL, Filter]

            # build up new INFO field
            # will not be affected by number of samples
            # pick out the allele'th value from each of the info values
            new_info_values = [x[allele] for x in INFO_values]
            # create empty list that will be the new info field
            new_info = [None] * (len(new_info_values) + len(INFO_names))
            # fill it with the info names and new values
            new_info[::2] = INFO_names
            new_info[1::2] = new_info_values
            # paste the name with the corresponding value (which will be directly
            # after the name in the list)
            si = iter(new_info)
            new_info = [x + next(si, '') for x in si]
            # print(new_info)
            # join all of the info fields with ;
            new_info = ";".join(new_info)

            # add new INFO field to new_line
            new_line.append(new_info)

            # add the FORMAT names to the new line
            new_FORMAT_name = ":".join(FORMAT_name)
            new_line.append(new_FORMAT_name)
            # print(new_line)
            # new_line = ['MT', '3103', '.', 'CTACN', 'CTAC', '68191.9', '.', 'DP=7575;MQM=70;MQMR=0;QA=75397;QR=0;SAF=1289;SAR=1351;SRF=0;SRR=0', 'GT:DP:AD:RO:QR:AO:QA']
            # new_line = ['MT', '3103', '.', 'CTACN', 'CTACT', '68191.9', '.', 'DP=7575;MQM=60;MQMR=0;QA=5623;QR=0;SAF=100;SAR=103;SRF=0;SRR=0', 'GT:DP:AD:RO:QR:AO:QA']

            # Add new FORMAT field
            # FORMAT field is affected by number of samples

            for samp in range(0, n_samples):
                temp_format = FORMAT_values[samp]
                logging.debug(temp_format)

                # split the allele values in the FORMAT
                temp_format = [x.split(",") for x in temp_format]
                # if temp_format[0] == ['2/2']:
                # print(temp_format)
                # print(temp_format)

                # GT will be overwritten later
                # DP will only have one value per sample, and will be unaffected by allele
                # RO will only have one value per sample, and will be unaffected by allele
                # QR will only have one value per sample, and will be unaffected by allele
                #
                # AD, AO and QA need to be updated - depend on allele

                # AO_idx is where AO sits in the FORMAT field
                AO_idx = FORMAT_name.index("AO")
                # pick out the allele'th value of AO
                new_AO = temp_format[AO_idx][allele]
                # print(new_AO)

                # update the format field
                # print(temp_format)
                temp_format[AO_idx] = [new_AO]
                # print(temp_format)
                # print("")

                # QA_idx is where QA sits in the FORMAT field
                QA_idx = FORMAT_name.index("QA")
                # pick out the allele'th value of QA
                new_QA = temp_format[QA_idx][allele]
                # update the format field
                temp_format[QA_idx] = [new_QA]

                # AD_idx is where AD sits in the FORMAT field
                AD_idx = FORMAT_name.index("AD")
                # pick out the allele'th + 1 value of AD
                # +1 because the first value is always the reference depth
                new_AD = temp_format[AD_idx][allele + 1]
                # add the reference depth in
                new_AD = [temp_format[AD_idx][0], new_AD]
                # print(new_AD)
                # update the format field
                temp_format[AD_idx] = new_AD
                # print(temp_format)
                temp_format = [",".join(x) for x in temp_format]
                # print(temp_format)
                temp_format = ":".join(temp_format)
                # print(temp_format)

                # add the new format field to the end of the new line
                new_line.append(temp_format)
                # print(new_line)

            split_vcf.append(new_line)

    return (split_vcf)


# variant_list is a list, with no multiallelic variants
# Eg variant_list = [['MT', '3103', '.', 'CTACN', 'CTAC', '68191.9', '.', 'DP=7575; ... ], ['MT', '3103', '.', 'CTACN', 'CTACT', '68191.9', '.', 'DP=7575;]]
def split_MNP(variant_list):
    table = []

    for line in variant_list:
        # print(line)
        ref = list(line[3])
        ref_length = len(ref)

        alt = list(line[4])
        alt_length = len(alt)

        if alt_length > 1 and ref_length > 1:
            # then this is a MNP

            if alt_length != ref_length:
                # then ref and alt dont have the ame number of bases
                # there isnt much we can do to simplify apart from trim the starting bases if they are matching
                # could be a complex variant if it doesnt simplify after trimming

                min_length = min(alt_length, ref_length)
                num_matching_bases = 0
                # only loop over the minimum of alt or ref bases
                # count how many of the bases at the start are the same
                for base in range(0, min_length):
                    if ref[base] != alt[base]:
                        break
                    else:
                        num_matching_bases += 1
                # trim_num is the number of bases at the start to trim
                trim_num = num_matching_bases - 1
                temp_line = list(line)
                if trim_num > 0:
                    # update the position
                    temp_line[1] = int(line[1]) + trim_num
                    # trim the reference
                    temp_line[3] = temp_line[3][trim_num:]
                    # trim the alternate
                    temp_line[4] = temp_line[4][trim_num:]

                # if the new reference and alternate length are still > 1 then
                # then this is defined as a complex variant
                if len(temp_line[3]) > 1 and len(temp_line[4]) > 1:
                    # print("its still a complex variant")
                    temp_line[7] = temp_line[7] + ";TYPE=complex"
                    # print(temp_line)

                table.append(temp_line)

            else:
                # the alternate and the reference have the same number of bases
                # so we just try to remove the extra bases
                # This logic could be changed so that variants next to each other are reported together
                # eg AT -> CG is currently reported A->C and T->G
                for base_idx in range(0, ref_length):

                    if alt[base_idx] != ref[base_idx]:
                        # then this is a variant
                        temp_line = list(line)
                        # update the position, reference and alternate
                        temp_line[1] = int(line[1]) + base_idx
                        temp_line[3] = ref[base_idx]
                        temp_line[4] = alt[base_idx]

                        table.append(temp_line)
        else:
            # If either the alt or reference is one base, then not an MNP
            table.append(line)

    return (table)


# variant_list is a list of sublists, each sublist is  line in the vcf
# variant_list is a list, with no multiallelic variants
# Eg variant_list = [['MT', 73.0, '.', 'A', 'G', '148165', '.', 'DP=5912...], ['MT', 146.0, '.', 'T', 'C', '198193', '.', 'DP=7697;...], ['MT', '182', '.', 'C', 'T', '6.00899e-14', '.', 'DP=7955;...]]
def combine_lines(variant_list, p=0.002):
    # make dictionary with keys chromosome, position and alternate

    # print(variant_list)
    # we are interested in the lines that have the same chromosome position and alternate -
    # to identify these lines, create an "id" for each line that is position_alternate
    # if any lines have matching ids then they should be combined.

    # get the chromosome positions and alternates

    chromo = [x[0] for x in variant_list]
    unique_chromosomes_temp = list(set(chromo))

    # trying to keep the variants in a rough order (set doesnt preserve order)
    autosome = [str(x) for x in list(range(1, 23))]
    unique_chromosomes = [x for x in autosome if x in unique_chromosomes_temp]
    unique_chromosomes = unique_chromosomes + [x for x in ['X', 'Y', 'MT'] if x in unique_chromosomes_temp]
    unique_chromosomes = unique_chromosomes + [x for x in unique_chromosomes_temp if x not in unique_chromosomes]

    # print(unique_chromosomes)
    # sys.exit()
    # unique_chromosomes = [str(x) for x in unique_chromosomes]
    # x = list(range(1,23))
    # print(x)
    new_vcf = []
    # we loop over the chromsomes to fill new_vcf
    for uniq_chrom in unique_chromosomes:
        logging.debug('Processing chromosome ' + str(uniq_chrom))
        chromo_variant_list = [line for line in variant_list if line[0] == uniq_chrom]
        # print(chromo_variant_list)
        # sys.exit()
        # TODO: could remove the lines from the variant list now - would make faster?
        # print(uniq_chrom)
        # print(chromo_variant_list)
        # sys.exit()
        position = [x[1] for x in chromo_variant_list]
        # print(position)
        alternate = [x[4] for x in chromo_variant_list]

        # paste the position and alternate to make id called pos_alt_id
        line_id = []
        for pos, alt in zip(position, alternate):
            line_id.append("_".join([str(pos), alt]))
        # print(line_id)
        # line_id = ['MT_73_G', 'MT_116_C', 'MT_139_G', 'MT_146_C', ... , 'MT_16263_C', 'MT_16264_T', 'MT_16263_C', 'MT_16270_T', ...]

        # line_number is the line number from chromo_variant_list
        # this is to keep track of whats already been added to the new vcf
        line_number = list(range(len(chromo_variant_list)))
        # print(line_number)
        # line_number = [0, 1, 2, 3, ... 497, 498, 499, 500]
        while len(line_id) > 0:
            # take the first id
            x = line_id[0]
            # print(x)
            # x = 73_G
            # or later in the loop
            # x = 16263_C

            # find all the matching ids, including the first
            # print("finding matching lines")
            matching_lines_idx = [i for i, j in enumerate(line_id) if j == x]
            # print("finished findinf matching lines")
            # print(matching_lines_idx)
            # matching_lines_idx = [0] i.e. there is only one 73_G
            # or furher into the loop
            # matching_lines_idx = [0,2] i.e. there are 2 16263_C

            matching_lines = [line_number[i] for i in matching_lines_idx]
            # print(matching_lines)
            # matching_lines = [0]
            # or later into the loop for 16263_C
            # matching_lines = [497, 499]

            # remove the line_id and the line_number that have been assigned to matching_lines
            for i in sorted(matching_lines_idx, reverse=True):
                del line_id[i]
            for i in sorted(matching_lines_idx, reverse=True):
                del line_number[i]
                # print(line_id)
            # line_id = ['116_C', '139_G', '146_C', '152_C', '182_T', .... ]
            # or later into the loop
            # if x == '16263_C':
            # print(line_id)
            # line_id = ['16264_T', '16270_T', '16274_A', '16278_T', '16284_G', '16288_C', '16290_T', '16293_C', '16301_T', '16311_C', '16343_C', ... ]
            # ie 16263_C has been removed everywhere in list

            # print(line_number)
            # line_number = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, ... ]
            # or later into the loop
            # if x == '16263_C':
            # print(line_number)
            # line_number = [498, 500, 501, 502, 503 ... ]

            # TODO: the way I manipulate the INFO line is different in each case
            # should probably make the same so its easier to follow
            if len(matching_lines) == 1:
                # print("not repeated position")
                # print(matching_lines)
                # This is a unique record for this genomic position.
                # we just need to update and add values to the info/format field.

                # pick out the line from chromo_variant_list:
                # need the [0] to remove it from list, better way to do this? ie [[x]] to [x]
                new_line = [chromo_variant_list[i] for i in matching_lines][0]
                # print(new_line)

                INFO = new_line[7]
                INFO = INFO.split(";")

                INFO_names = [i.split("=")[0] for i in INFO]
                # print(INFO_names)
                INFO_values = [i.split("=")[1] for i in INFO]

                # we want to caclulate SBR and add to INFO
                # SBR = strand bias of reference
                # SBR = SRF/(SRF+SRR)

                SRF_idx = INFO_names.index('SRF')
                SRF = float(INFO_values[SRF_idx])

                SRR_idx = INFO_names.index('SRR')
                SRR = float(INFO_values[SRR_idx])

                # to avoid dividing by zero:
                if SRF > 0 or SRR > 0:
                    SBR = round(SRF / (SRF + SRR), 3)
                else:
                    SBR = 0

                INFO_values.append(SBR)
                INFO_names.append("SBR")

                # we want to calculate SBA and add to INFO
                # SBA = strand bias of alternate
                # i.e. SBA = SAF/(SAF+SAR)

                SAF_idx = INFO_names.index('SAF')
                SAF = float(INFO_values[SAF_idx])

                SAR_idx = INFO_names.index('SAR')
                SAR = float(INFO_values[SAR_idx])
                # note SAF will always > 0, otherwise variant wouldn't be reported
                SBA = round(SAF / (SAF + SAR), 3)

                INFO_values.append(SBA)
                INFO_names.append("SBA")

                # we want to add VAF to the FILTER field
                # to do this we need AO and DP
                # we also want to update the genotype
                # and add AQR - average base quality of the reference
                # we update AD to be RO,AO, because of the bug in earlier version that had AD wrong, just to be safe
                # these are all per sample

                # because not matching position, there will only be one line of FORMAT values
                FORMAT_values = new_line[9:]
                FORMAT_values = [x.split(":") for x in FORMAT_values]
                # print(FORMAT_values)
                n_samples = len(FORMAT_values)

                FORMAT_names = new_line[8]
                FORMAT_names = FORMAT_names.split(":")

                # add AQR after QR
                QR_idx = FORMAT_names.index('QR')
                FORMAT_names.insert((QR_idx + 1), 'AQR')
                # add AQA after QA
                QA_idx = FORMAT_names.index('QA')
                FORMAT_names.insert((QA_idx + 1), 'AQA')
                # print(FORMAT_names)

                # add VAF to the format
                FORMAT_names.append('VAF')
                # add q to the formal
                FORMAT_names.append('q')
                # print(FORMAT_names)
                # add dummy VAF placeholder
                [x.extend(['-1']) for x in FORMAT_values]
                # add dummy q placeholder
                [x.extend(['-1']) for x in FORMAT_values]
                # print(FORMAT_values)

                AO_vector = []
                # VAF = []
                DP_vector = []
                # new_genotype = []
                for samp in range(0, n_samples):
                    # print("")
                    temp_format = FORMAT_values[samp]
                    # print(temp_format)
                    # add dummy placeholder for AQR and AQA
                    temp_format.insert((QR_idx + 1), '-1')
                    temp_format.insert((QA_idx + 1), '-1')
                    # print(temp_format)

                    AO_idx = FORMAT_names.index('AO')
                    AO = int(temp_format[AO_idx])
                    AO_vector.append(AO)
                    # print(AO)
                    # print(AO_vector)

                    RO_idx = FORMAT_names.index('RO')
                    RO = int(temp_format[RO_idx])

                    # change the format DP to be AO + RO
                    DP_idx = FORMAT_names.index('DP')
                    DP = RO + AO
                    temp_format[DP_idx] = str(DP)
                    DP_vector.append(DP)
                    # DP = int(temp_format[DP_idx])
                    # print(DP)

                    # DP = AO + RO can be zero if another sample cal
                    if DP != 0:
                        samp_VAF = round(float(AO) / float(DP), 4)
                    else:
                        samp_VAF = 0

                    # VAF.append(str(samp_VAF))
                    VAF_idx = FORMAT_names.index('VAF')
                    temp_format[VAF_idx] = str(samp_VAF)

                    # add in q
                    q = mity_qual(AO, DP, p=p)
                    q_idx = FORMAT_names.index('q')
                    temp_format[q_idx] = str(q)

                    # print(p)
                    # print(float(DP))
                    # print(float(AO))
                    # print(q)
                    # exit()
                    # add in the new genotype
                    if AO < 4:
                        new_genotype = "0/0"
                    else:
                        if samp_VAF == 1:
                            new_genotype = "1/1"
                        elif samp_VAF < 1:
                            new_genotype = "0/1"

                    # print(VAF)
                    # if samp_VAF > 0.9:
                    #   new_genotype = "1/1"
                    # elif samp_VAF > 0:
                    #   new_genotype = "0/1"
                    # else:
                    #   new_genotype = "0/0"
                    # print(new_genotype)

                    # replace the genotypes
                    GT_idx = FORMAT_names.index('GT')
                    temp_format[GT_idx] = new_genotype
                    # print(temp_format)
                    # print("")

                    RO_idx = FORMAT_names.index('RO')
                    RO = int(temp_format[RO_idx])
                    # print(RO)

                    # calculate new AD
                    new_AD = [str(RO), str(AO)]
                    new_AD = ",".join(new_AD)
                    AD_idx = FORMAT_names.index('AD')
                    temp_format[AD_idx] = new_AD
                    # print(temp_format)

                    # calculate AQR
                    QR_idx = FORMAT_names.index('QR')
                    # print(AO_idx)
                    QR = float(temp_format[QR_idx])
                    if RO > 0:
                        AQR = str(round(QR / RO, 3))
                    else:
                        # if there are no reference reads, set AQR high so you dont filter on it
                        # AQR = 1000000
                        # AQR = 1000000 will be confusing, so check in the filters function if RO = 0
                        AQR = str(0)
                    # add AQR to the FORMAT field
                    AQR_idx = FORMAT_names.index('AQR')
                    temp_format[AQR_idx] = AQR

                    # calculate AQA
                    QA_idx = FORMAT_names.index('QA')
                    QA_idx = FORMAT_names.index('QA')
                    # print(AO_idx)
                    QA = float(temp_format[QA_idx])
                    # here AO could be less than zero because it could be another sample that has this variant
                    if AO > 0:
                        AQA = str(round(QA / AO, 3))
                    else:
                        AQA = str(0)
                    # add AQA to the FORMAT field
                    AQA_idx = FORMAT_names.index('AQA')
                    temp_format[AQA_idx] = str(AQA)

                    # add the updated FORMAT field for the given sample to the new_line
                    temp_format = ":".join(temp_format)
                    new_line[samp + 9] = temp_format

                # update QUAL
                cohort_AO = sum(int(x) for x in AO_vector)
                cohort_DP = sum(int(x) for x in DP_vector)
                cohort_QUAL = mity_qual(cohort_AO, cohort_DP, p=p)
                new_line[5] = str(round(cohort_QUAL, 1))
                logging.debug("Cohort QUAL score: {} := {} vs {}".format(new_line[5], ",".join(str(x) for x in AO_vector), ",".join(str(x) for x in DP_vector)))

                # add the new FORMAT names
                FORMAT_names = ":".join(FORMAT_names)
                new_line[8] = FORMAT_names
                # print(new_line)

                # add VAF to the info
                # VAF = ",".join(VAF)
                # INFO_values.append(VAF)
                # INFO_names.append("VAF")

                # dont update the INFO depth - doesnt make sense
                # when you have more than one sample
                # unless you do RO+AO for each sample
                # ie total depth at that position in all the samples...

                INFO_DP_idx = INFO_names.index('DP')
                new_info_DP = sum(DP_vector)
                INFO_values[INFO_DP_idx] = new_info_DP
                # print(DP_vector)
                # print(new_info_DP)

                # put the = sign back
                INFO_names = [x + "=" for x in INFO_names]
                new_info = [None] * (len(INFO_names) + len(INFO_values))
                new_info[::2] = INFO_names
                new_info[1::2] = INFO_values
                new_info = [str(x) for x in new_info]
                si = iter(new_info)
                new_info = [x + next(si, '') for x in si]
                new_info = ";".join(new_info)
                # print(new_info)
                # update INFO
                new_line[7] = new_info
                # print(new_line)
                new_vcf.append(new_line)


            else:
                # then this is a repeated position
                # we need to update and add values to the info/format field, as well as combine the matching lines
                # print("repeated position")
                # print(matching_lines)

                new_line = [chromo_variant_list[i] for i in matching_lines]
                replacement_line = new_line[0]
                # print(new_line)

                #######################
                ##### FORMAT field
                new_FORMAT = []

                # loop over samples, and update each sample separately and then add to new format
                # each line should have same number of samples
                # TODO could move the calculation of the number of samples to the start of the function
                first_FORMAT_values = new_line[0][9:]
                n_samples = len(first_FORMAT_values)

                # FORMAT_values = new_line[9:]
                # FORMAT_values = [x.split(":") for x in FORMAT_values]
                # print(FORMAT_values)

                FORMAT_names = new_line[0][8]
                FORMAT_names = FORMAT_names.split(":")
                # print(FORMAT_names)

                # we need the AO vector later for weighted mean
                AO_idx = FORMAT_names.index('AO')
                AO_sum_vector = []
                for line in new_line:
                    # print(line)
                    formats = line[9:]
                    formats = [x.split(":") for x in formats]
                    AOs = [float(x[AO_idx]) for x in formats]
                    AO_sum_vector.append(sum(AOs))
                # print(AO_sum_vector)

                # add AQR after QR
                QR_idx = FORMAT_names.index('QR')
                FORMAT_names.insert((QR_idx + 1), 'AQR')
                # add AQA after QA
                QA_idx = FORMAT_names.index('QA')
                FORMAT_names.insert((QA_idx + 1), 'AQA')
                # print(FORMAT_names)
                # add the updated FORMAT_names to the new_line

                # add VAF to the format
                FORMAT_names.append('VAF')
                # add q to the format
                FORMAT_names.append('q')
                # print(FORMAT_names)
                # add dummy VAF placeholder
                # [x.extend(['-1']) for x in FORMAT_values]
                # print(FORMAT_values)

                new_AD = []
                VAF = []
                # because we need it later, create and keep an new_AO vector
                new_AO_vector = []
                new_DP_vector = []
                for samp in range(0, n_samples):
                    #   # print("")
                    temp_format = [i[9 + samp].split(":") for i in new_line]
                    # print(temp_format)

                    # add dummy placeholder for AQR and AQA

                    for x in range(0, len(temp_format)):
                        temp_format[x].insert((QR_idx + 1), '-1')
                        temp_format[x].insert((QA_idx + 1), '-1')
                    # print(temp_format)

                    # AO
                    AO_idx = FORMAT_names.index('AO')
                    # TODO: do we need new_AO_vector anymore?
                    AO_vector = [i[AO_idx] for i in temp_format]
                    new_AO = str(int(Sum(AO_vector)))
                    new_AO_vector.append(new_AO)
                    # new_AO = str(int(new_AO))
                    # print(AO_vector)
                    # print(new_AO)

                    # RO
                    RO_idx = FORMAT_names.index('RO')
                    RO_vector = [i[RO_idx] for i in temp_format]
                    new_RO = unchanged(RO_vector)
                    new_RO = str(int(new_RO))
                    # print(RO_vector)
                    # print(new_RO)

                    new_AD = [str(new_RO), str(new_AO)]
                    new_AD = ",".join(new_AD)

                    # DP
                    DP_idx = FORMAT_names.index('DP')
                    new_DP = str(int(new_AO) + int(new_RO))
                    new_DP_vector.append(new_DP)
                    # DP_vector = [i[DP_idx] for i in temp_format]
                    # new_DP = unchanged(DP_vector)
                    # new_DP = str(int(new_DP))
                    # print(DP_vector)
                    # print(new_DP)

                    if new_DP != 0:
                        samp_VAF = round(float(new_AO) / float(new_DP), 4)
                    else:
                        samp_VAF = 0.0
                    logging.debug("Variant has samp_VAF: {}".format(str(samp_VAF)))
                    VAF.append(str(samp_VAF))
                    # print(VAF)

                    samp_q = mity_qual(new_AO, new_DP, p=p)
                    samp_q = str(samp_q)

                    # update GT
                    if int(new_AO) < 4:
                        new_genotype = "0/0"
                    else:
                        if samp_VAF >= 0.95:
                            new_genotype = "1/1"
                        elif samp_VAF < 0.95:
                            new_genotype = "0/1"
                    samp_VAF = str(samp_VAF)

                    # QR
                    QR_vector = [i[FORMAT_names.index('QR')] for i in temp_format]
                    new_QR = unchanged(QR_vector)
                    new_QR = str(new_QR)

                    # calculate AQR
                    if int(new_RO) > 0:
                        AQR = str(float(new_QR) / float(new_RO))
                    else:
                        # if there are no reference reads, set AQR high so you dont filter on it
                        # AQR = 1000000
                        # AQR = 1000000 will be confusing, so check in the filters function if RO = 0
                        AQR = str(0)

                    # QA
                    QA_idx = FORMAT_names.index('QA')
                    QA_vector = [i[QA_idx] for i in temp_format]
                    new_QA = Sum(QA_vector)
                    new_QA = str(int(new_QA))

                    # calculate AQA
                    # here AO could be less than zero because it could be another sample that has this variant
                    if int(new_AO) > 0:
                        AQA = str(float(new_QA) / float(new_AO))
                    else:
                        # if there are no reference reads, set AQR high so you dont filter on it
                        # AQR = 1000000
                        # AQR = 1000000 will be confusing, so check in the filters function if RO = 0
                        AQA = str(0)

                    new_FORMAT.append(":".join(
                        [new_genotype, new_DP, new_AD, new_RO, new_QR, AQR, new_AO, new_QA, AQR, samp_VAF, samp_q]))
                    # print(new_FORMAT)

                # print(new_FORMAT)

                #######################
                ##### INFO field
                INFO_field = [i[7].split(":") for i in new_line]
                # print(INFO_field)
                INFO_field = [x[0].split(";") for x in INFO_field]
                # print(INFO_field)

                INFO_values = []
                for line in INFO_field:
                    INFO_values.append([x.split("=")[1] for x in line])

                # print(INFO_values)

                INFO_names = INFO_field[0]
                INFO_names = [x.split("=")[0] for x in INFO_names]
                # print(INFO_names)

                # DP
                DP_idx = INFO_names.index('DP')
                # print(DP_idx)
                DP_vector = [i[DP_idx] for i in INFO_values]
                # print(DP_vector)
                new_DP = unchanged(DP_vector)
                # new_DP = sum(DP_vector)
                new_DP = "DP=" + str(new_DP)
                # print(new_DP)

                # MQM - Mean mapping quality of observed alternate alleles
                MQM_idx = INFO_names.index('MQM')
                # print(MQM_idx)
                MQM_vector = [i[MQM_idx] for i in INFO_values]
                new_MQM = weighted_mean(mean_list=MQM_vector, count_list=AO_sum_vector)
                new_MQM = "MQM=" + str(new_MQM)
                # print(new_MQM)

                # MQMR - Mean mapping quality of observed reference alleles
                MQMR_idx = INFO_names.index('MQMR')
                MQMR_vector = [i[MQMR_idx] for i in INFO_values]
                # print(MQMR_vector)
                new_MQMR = unchanged(MQMR_vector)
                new_MQMR = "MQMR=" + str(new_MQMR)
                # print(new_MQMR)

                # QA - Sum of quality of the alternate observations
                QA_idx = INFO_names.index('QA')
                QA_vector = [i[QA_idx] for i in INFO_values]
                new_QA = Sum(QA_vector)
                new_QA = "QA=" + str(new_QA)

                # QR - Sum of quality of the reference observations
                QR_idx = INFO_names.index('QR')
                QR_vector = [i[QR_idx] for i in INFO_values]
                new_QR = unchanged(QR_vector)
                new_QR = "QR=" + str(new_QR)

                # SAF - Number of alternate observations on the forward strand
                SAF_idx = INFO_names.index('SAF')
                SAF_vector = [i[SAF_idx] for i in INFO_values]
                new_SAF = Sum(SAF_vector)

                # SAR - Number of alternate observations on the reverse strand
                SAR_idx = INFO_names.index('SAR')
                SAR_vector = [i[SAR_idx] for i in INFO_values]
                new_SAR = Sum(SAR_vector)

                SBA = new_SAR / (new_SAR + new_SAF)

                # SRF - Number of reference observations on the forward strand
                SRF_idx = INFO_names.index('SRF')
                SRF_vector = [i[SRF_idx] for i in INFO_values]
                new_SRF = unchanged(SRF_vector)
                # print(new_SRF)

                # SRR - Number of reference observations on the forward strand
                SRR_idx = INFO_names.index('SRR')
                SRR_vector = [i[SRR_idx] for i in INFO_values]
                new_SRR = unchanged(SRR_vector)
                #  TODO: why do i need int here?
                if int(new_SRF) > 0 or int(new_SRR) > 0:
                    SBR = int(new_SRF) / (int(new_SRF) + int(new_SRR))
                else:
                    SBR = 0

                new_SAF = "SAF=" + str(new_SAF)
                new_SAR = "SAR=" + str(new_SAR)
                new_SRF = "SRF=" + str(new_SRF)
                new_SRR = "SRR=" + str(new_SRR)
                new_SBA = "SBA=" + str(SBA)
                new_SBR = "SBR=" + str(SBR)

                if "TYPE" in INFO_names:
                    # then add normalised to it
                    new_type = "TYPE=complex,normalised"
                else:
                    new_type = "TYPE=normalised"

                new_info = [new_DP, new_MQM, new_MQMR, new_QA, new_QR, new_SAF, new_SAR, new_SRF, new_SRR, new_type,
                            new_SBA, new_SBR]
                # add the VAF to the end of the info line
                # VAF = ",".join(VAF)
                # VAF = "VAF="+str(VAF)
                # new_info.append(VAF)
                new_info = ";".join(new_info)

                cohort_AO = sum(int(x) for x in new_AO_vector)
                cohort_DP = sum(int(x) for x in new_DP_vector)
                cohort_QUAL = mity_qual(cohort_AO, cohort_DP, p=p)
                logging.debug("Cohort QUAL score: {} := {} vs {}".format(replacement_line[5], ",".join(str(x) for x in new_AO_vector), ",".join(str(x) for x in new_DP_vector)))
                replacement_line[5] = str(round(cohort_QUAL, 1))
                replacement_line[7] = new_info
                replacement_line[9:] = new_FORMAT
                FORMAT_names = ":".join(FORMAT_names)
                replacement_line[8] = FORMAT_names
                # print(replacement_line)
                new_vcf.append(replacement_line)
    # print(new_vcf)
    return (new_vcf)


def mity_qual(AO, DP, p=0.002):
    """
    Compute variant quality
    :param AO: (int) number of alternative reads
    :param DP: (int) total read depth
    :param p: (float) noise parameter. default=0.002
    :return: (float) phred-scaled quality score

    >>> [mity_qual(x,10) for x in [1,2,3,4,5,6,7,8,9,10]]
    '[37.49, 60.22, 84.78, 110.97, 138.75, 168.16, 199.4, 232.92, 269.9, 296.89]'
    >>> [mity_qual(x,100) for x in [5,10,15,20,25,30,35,40,45,50]]
    '[71.87, 156.08, 251.23, 354.34, 463.9, 579.05, 699.21, 824.04, 953.32, 1086.94]'
    >>> [mity_qual(x,1000) for x in [5,10,15,20,25,30,35,40,45,50]]
    '[17.84, 50.97, 93.59, 142.89, 197.36, 256.02, 318.25, 383.57, 451.61, 522.1]'
    >>> [mity_qual(x,10000) for x in [5,10,15,20,25,30,35,40,45,50]]
    '[0.0, 0.05, 0.74, 3.56, 9.51, 18.73, 31.0, 46.04, 63.57, 83.37]'
    >>> mity_qual(118,118)
    '3211.76'
    >>> mity_qual(119,119)
    '3220'
    >>> mity_qual(10000,10000)
    '3220'
    """
    q = 0.0
    AO = int(AO)
    DP = int(DP)
    if AO > 0 and DP > 0:
        if DP == AO:
            # then penalise homoplasmic variants with low numbers of reads
            DP = DP + 1
        # v1: the cdf implementation caps at 159.55, due to the cdf capping at 0.9999999999999999
        # q = round(abs(-10 * log10(1 - binom.cdf(AO, DP, p))), 2)
        # v2: the logsf implementation is identical to cdf for low values and continues to scale up
        with np.errstate(divide='ignore'):
            q = float(round(-4.342945 * binom.logsf(AO, DP, p), 2))
    if isinf(q):
        # this is the maximum q that we observed looking at homoplastic variants around ~129x depth
        q = float(3220)
    return q

def add_filter(variant_list, min_DP=15, SB_range=[0.1, 0.9], min_MQMR=30, min_AQR=20):
    """
    FILTER out poor quality variants in a VCF.
    If one samples in the vcf passes, then they all pass in the FILTER column.

    :param variant_list: the parsed VCF
    :param min_DP: minimum total read depth of the variant
    :param SB_range: min and max SB range, where SB is the strand bias of either the REF (SBR_FAIL), or ALT (SBA_FAIL) reads.
    :param min_MQMR: minimum MQMR, where MQMR is the "Mean mapping quality of observed reference alleles"
    :param min_AQR: minimum MQMR, where MQMR is the "Mean mapping quality of observed reference alleles"
    :return: a
    """
    # @TODO: refactor this to use pyvcf

    # print(variant_list)
    # sys.exit()

    #############################################
    ##### HARD CODED FILTER VALUES
    blacklist = list(range(302, 319))
    blacklist = blacklist + [3105, 3106, 3107, 3108]

    new_vcf = []

    for line in variant_list:
        # print(line)

        POS = int(line[1])

        info = line[7].split(";")

        info_names = [x.split("=")[0] for x in info]
        info_values = [x.split("=")[1] for x in info]

        SBR = float(info_values[info_names.index("SBR")])
        SBA = float(info_values[info_names.index("SBA")])
        MQMR = float(info_values[info_names.index("MQMR")])

        FORMAT_names = line[8].split(":")

        FORMAT = line[9:]

        samp_POS_fil = []
        samp_SBR_fil = []
        samp_SBA_fil = []
        samp_MQMR_fil = []
        samp_AQR_fil = []
        # @clare, why don't you just fail variants with depth < min_DP in one filter (samp_DP_fil)?

        for samp in FORMAT:
            samp_FORMAT = samp.split(":")

            AO = int(samp_FORMAT[FORMAT_names.index("AO")])
            RO = int(samp_FORMAT[FORMAT_names.index("RO")])
            AQR = float(samp_FORMAT[FORMAT_names.index("AQR")])

            from .util import make_hgvs
            hgvs = make_hgvs(int(line[1]), line[3], line[4])
            if RO > min_DP:
                # logging.debug("Position {} pass min_DP filter, with RO {} > {}".format(hgvs, str(RO), str(min_DP)))
                if not SB_range[0] <= SBR <= SB_range[1]:
                    logging.debug("Position {} failed SBR filter {}".format(hgvs, str(SBR)))
                    samp_SBR_fil.append(0)

                if MQMR < min_MQMR:
                    samp_MQMR_fil.append(0)

                if AQR < min_AQR:
                    samp_AQR_fil.append(0)

            if AO > min_DP:
                # logging.debug("Position {} pass min_DP filter, with AO {} > {}".format(hgvs, str(AO), str(min_DP)))
                if not SB_range[0] <= SBA <= SB_range[1]:
                    logging.debug("Position {} failed SBA filter {}".format(hgvs, str(round(SBA, 3))))
                    samp_SBA_fil.append(0)

            if POS in blacklist:
                samp_POS_fil.append(0)

        # if they aren't all zero than at least one passed
        # only one sample has to pass
        FILTER = []
        if len(samp_SBR_fil) < len(FORMAT):
            SBR_fil = 1
        else:
            SBR_fil = 0
            FILTER.append("SBR")

        if len(samp_SBA_fil) < len(FORMAT):
            SBA_fil = 1
        else:
            SBA_fil = 0
            FILTER.append("SBA")

        if len(samp_MQMR_fil) < len(FORMAT):
            MQMR_fil = 1
        else:
            MQMR_fil = 0
            FILTER.append("MQMR")

        if len(samp_AQR_fil) < len(FORMAT):
            AQR_fil = 1
        else:
            AQR_fil = 0
            FILTER.append("AQR")

        if len(samp_POS_fil) < len(FORMAT):
            POS_fil = 1
        else:
            POS_fil = 0
            FILTER.append("POS")

        # FILTER = "FAIL"
        if sum([SBR_fil, SBA_fil, MQMR_fil, AQR_fil, POS_fil]) == 5:
            FILTER = "PASS"
        else:
            FILTER = ";".join(FILTER)
        # print(FILTER)

        line[6] = FILTER

        # line.append(FILTER_INFO)
        # print(line)
        # if POS == 302:
        #   print(info)
        #   print(POS_fil)
        #   print(line)
        new_vcf.append(line)
    # print(new_vcf)
    return (new_vcf)


def update_header(col_names, header_lines, p, SB_range=[0.1, 0.9], min_MQMR=30, min_AQR=20):
    ##############################
    ####### Process Header lines
    header_lines.append([
        '##ALT=<ID=NON_REF,Description="Represents any '
        'possible alternative allele at this location">'])
    ##############################
    ######## INFO lines
    # DP - across all the samples. Is this then summed across all the samples?
    header_lines.append([
        '##INFO=<ID=DP,Number=1,Type=Integer,'
        'Description="Approximate read depth; some reads may '
        'have been filtered">'])
    # MQM - according to freebayes there should be one per allele (Number=A).
    header_lines.append([
        '##INFO=<ID=MQM,Number=A,Type=Float,Description="Mean '
        'mapping quality of observed alternate alleles">'])
    # MRMR only ever one. Could also add to format field, see MQM
    header_lines.append([
        '##INFO=<ID=MQMR,Number=1,Type=Float,'
        'Description="Mean mapping quality of observed '
        'reference alleles">'])
    # QA Again one per allele. But there are sample specific alleles in the
    # format field
    header_lines.append([
        '##INFO=<ID=QA,Number=A,Type=Integer,'
        'Description="Alternate allele quality sum in phred">'])
    # QR only one per reference. Sample specific in the format field
    header_lines.append([
        '##INFO=<ID=QR,Number=1,Type=Integer,'
        'Description="Reference allele quality sum in phred">'])
    # SAF, SAR, SRF, SRR - onle per allele, averaged over all samples. Consider
    # adding sample specific SBR and SBA
    header_lines.append([
        '##INFO=<ID=SAF,Number=A,Type=Integer,'
        'Description="Total number of alternate observations '
        'on the forward strand">'])
    header_lines.append([
        '##INFO=<ID=SAR,Number=A,Type=Integer,'
        'Description="Total number of alternate observations '
        'on the reverse strand">'])
    header_lines.append([
        '##INFO=<ID=SBA,Number=1,Type=Float,'
        'Description="Strand bias of the alternate reads, '
        'SBA=SAF/(SAF+SAR)">'])
    header_lines.append([
        '##INFO=<ID=SBR,Number=1,Type=Float,'
        'Description="Strand bias of the reference reads, '
        'SBR=SRF/(SRF+SRR)">'])
    header_lines.append([
        '##INFO=<ID=SRF,Number=1,Type=Integer,'
        'Description="Number of reference observations on the '
        'forward strand">'])
    header_lines.append([
        '##INFO=<ID=SRR,Number=1,Type=Integer,'
        'Description="Number of reference observations on the '
        'reverse strand">'])
    header_lines.append([
        '##INFO=<ID=TYPE,Number=1,Type=String,'
        'Description="normalised := when the variant has been normalised, complex := when the new reference and alternate length are still > 1 after normalising">'])
    # vaf - will be one for each allele. it is very useful so its good to be in
    # the INFO field
    ##############################
    ######## FILTER lines:
    header_lines.append([
        '##FILTER=<ID=POS,Description="Variant falls in '
        'the blacklist of positions: MT:302-319, '
        'MT:3105-3108">'])
    header_lines.append([
        '##FILTER=<ID=SBR,Description="For all alleles RO '
        '> 15 and (SBR > ' + str(SB_range[1]) + ' or SBR < ' + str(SB_range[0]) + ')">'])
    header_lines.append([
        '##FILTER=<ID=SBA,Description="For all alleles AO '
        '> 15 and (SBA > ' + str(SB_range[1]) + ' or SBA < ' + str(SB_range[0]) + ')">'])
    header_lines.append(
        ['##FILTER=<ID=MQMR,Description="For all alleles MQMR<' + str(min_MQMR) + '">'])
    header_lines.append(
        ['##FILTER=<ID=AQR,Description="For all alleles AQR<' + str(min_AQR) + '">'])
    # header_lines.append(['##FILTER=<ID=VAF,Number=A,Type=Float,
    # Description="Allele frequency in the range (0,1] - the ratio of the
    # number of alternate reads to reference reads">'])
    ##############################
    ######## FORMAT lines
    # GT, DP, RO, AO, QR, QA, AD
    header_lines.append([
        '##FORMAT=<ID=AD,Number=.,Type=Integer,'
        'Description="Allelic depths for the ref and alt '
        'alleles in the order listed">'])
    ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the
    # ref and alt alleles in the order listed">
    header_lines.append([
        '##FORMAT=<ID=AO,Number=1,Type=String,'
        'Description="Alternate observations">'])
    header_lines.append([
        '##FORMAT=<ID=AQA,Number=A,Type=Float,'
        'Description="Average base quality of the alternate '
        'reads, AQA=QA/AO">'])
    header_lines.append([
        '##FORMAT=<ID=AQR,Number=A,Type=Float,'
        'Description="Average base quality of the reference '
        'reads, AQR=QR/RO">'])
    header_lines.append([
        '##FORMAT=<ID=DP,Number=1,Type=Integer,'
        'Description="Read depth, AO+RO">'])
    header_lines.append([
        '##FORMAT=<ID=GT,Number=1,Type=String,'
        'Description="Genotype. Variants with less than 4 alt reads are 0/0, '
        'then variants with VAF < 0.95 are 0/1, and variants with VAF >= 0.95 are 1/1">'])
    header_lines.append([
        '##FORMAT=<ID=QA,Number=A,Type=Integer,Description="Alternate allele quality sum in phred">'])
    header_lines.append([
        '##FORMAT=<ID=QR,Number=1,Type=Integer,Description="Reference allele quality sum in phred">'])
    header_lines.append([
        '##FORMAT=<ID=RO,Number=1,Type=String,Description="Reference observations">'])
    header_lines.append([
        '##FORMAT=<ID=VAF,Number=A,Type=Float,Description="Allele frequency in the range (0,1] - the ratio of the number of alternate reads to reference reads">'])
    header_lines.append([
        '##FORMAT=<ID=q,Number=A,Type=Float,Description="Phred scaled binomial probability of seeing AO reads from DP, assuming a noise floor of p=' + str(
            p) + '. ">'])

    ##############################
    ######## Add col names back in
    header_lines.append([col_names])


def do_normalise(vcf, out_file=None, p=0.002, SB_range=[0.1,0.9], min_MQMR=30, min_AQR=20, chromosome=None, genome="mitylib/reference/hs37d5.genome"):
    """
    Normalise and FILTER a mity VCF file.

    This function splits multi-allelic variants and MNPs, whilst taking care to
    split the VCF header properly. After normalising, it also updates the FILTER,
    to remove poor quality variants, and calculates QUAL and GENOTYPE:QUAL scores.

    :param vcf: the path to a vcf.gz file. It can be single-sample,
    or multi-sample, and
         should have been generated by mity call, or FreeBayes
    :type vcf: str

    :param out_file: The name of the normalised vcf.gz file
    :type out_file: str

    :param p: The minimum noise threshold.
    :type p: float

    :param chromosome: discard variants not on this chromosome. If None,
    then all variants will be processed.
    :type chromosome: str

    :returns: Nothing. This creates a vcf.gz named out_file
    :rtype: None
    """
    if out_file is None:
        out_file = vcf.replace(".vcf.gz", ".norm.vcf.gz")

    file = gzip.open(vcf, 'rt')

    # split the header and the variants into two seperate lists
    # @TODO: refactor this to use pyvcf
    logging.debug('Splitting header and variants')
    col_names = None
    header_lines = []
    variants = []
    for line in file:
        if line[0] == "#":
            if re.match("#CHROM", line):
                col_names = line.strip()
                continue
            if re.match("##FORMAT", line) or re.match("##INFO", line):
                continue
            # only keep the non-FORMAT and non-INFO lines. We'll add these back later.
            line = line.strip()
            header_lines.append([line])
        else:
            line = line.strip()
            line_chromosome = line.split('\t')[0]
            if chromosome is None or line_chromosome == chromosome:
                variants.append(line)

    if len(variants) == 0:
        logging.warning('No variants in VCF with specified chromosome/s')

    logging.debug('Splitting multiallelic')
    single_allele = split_multi_allelic(variants)

    logging.debug('Splitting MNPs')
    no_mnp = split_MNP(single_allele)
    # debug_print_vcf_lines(no_mnp)
    logging.debug('Combining duplicated variants')
    combined_variants = combine_lines(no_mnp, p=p)
    # debug_print_vcf_lines(combined_variants)
    logging.debug('Adding filter variants')
    sba_filter = []
    filtered_variants = add_filter(combined_variants, SB_range=SB_range, min_MQMR=min_MQMR, min_AQR=min_AQR)
    # debug_print_vcf_lines(filtered_variants)
    update_header(col_names, header_lines, p=p, SB_range=SB_range, min_MQMR=min_MQMR, min_AQR=min_AQR)
    # debug_print_vcf_lines(filtered_variants)
    logging.info('Writing normalised vcf: {}'.format(out_file))

    new_vcf = header_lines + filtered_variants

    write_vcf(new_vcf, out_file, genome)


def debug_print_vcf_lines(x):
    for line in x:
        print(*line, sep='\t')
    sys.exit()
