import sys
import logging
import gzip
import pandas
import os.path
import xlsxwriter
from .util import check_missing_file, create_prefix, make_hgvs, get_annot_file

def make_table(variants, samples, vep_headers, impact_dict, min_vaf):
    table = []
    n_sample = len(samples)
    
    for lineparts in variants:
        # These should be the same for all VCFs
        chromosome = lineparts[0]
        pos = lineparts[1]
        ref = lineparts[3]
        alt = lineparts[4]
        QUAL = lineparts[5]
        FILTER = lineparts[6]

        hgvs = make_hgvs(pos, ref, alt)
        # print(hgvs)
        info = lineparts[7]
        info = info.split(";")
        # print(info)
        
        ## get the filters
        
        POS_FILTER = "FAIL" if 'POS' in FILTER else "PASS"
        SBR_FILTER = "FAIL" if 'SBR' in FILTER else "PASS"
        SBA_FILTER = "FAIL" if 'SBA' in FILTER else "PASS"
        MQMR_FILTER = "FAIL" if 'MQMR' in FILTER else "PASS"
        AQR_FILTER = "FAIL" if 'AQR' in FILTER else "PASS"

        INFO_names = [x.split("=")[0] for x in info]
        INFO_values = [i.split("=")[1] for i in info]
        # print(INFO_names)
        # print(INFO_values)
        
        ############################################
        ########## VEP
        if vep_headers != "":
            
            VEP_idx = INFO_names.index('CSQ')
            VEP = INFO_values[VEP_idx]
            VEP = VEP.split(",")
            VEP = [v.split('|') for v in VEP]
            # print(VEP)
            
            # get the index of where the consequence sits in the list of VEP 
            # annotations
            consequence_idx = vep_headers.index('Consequence')
            
            # get the index of the annotations that don't contain the word stream
            not_stream = ['stream' not in x[consequence_idx] for x in VEP]
            stream_idx = [i for i, x in enumerate(not_stream) if x]
            # print(stream_idx)
            
            # remove the sublists that have 'stream'in the consequence from 
            # the VEP annotations
            VEP = [VEP[i] for i in stream_idx]
            # print(VEP)
            # remove empty sublists
            VEP = [x for x in VEP if x != []]
            # print(VEP)
            
            # if there are no annotations left, completely remove the sublist
            if VEP == []:
                # fill with dot
                VEP = [['.'] * len(vep_headers)]
            else:
                VEP = [[element or '.' for element in sublist] for sublist in
                       VEP]
            
            # print(VEP)
            
            # get all the consequences for each variant in a single list
            # so that we can see what the highest impact is
            csq = [x[consequence_idx] for x in VEP]
            # some of them are split by '&'
            csq = [x.split("&") for x in csq]
            # unlist
            csq = sum(csq, [])
            # print(csq)
            
            # print(impact)
            if csq == ['.']:
                # impact = ["."]
                highest_impact = "."
            else:
                impact = [impact_dict[x] for x in csq]
                if 'HIGH' in impact:
                    highest_impact = 'HIGH'
                elif 'MODERATE' in impact:
                    highest_impact = 'MODERATE'
                elif 'MODIFIER' in impact:
                    highest_impact = 'MODIFIER'
                elif 'LOW' in impact:
                    highest_impact = 'LOW'
                else:
                    sys.error("unknown SO term")
                    
                    # print(highest_impact)
            
            # concatenate the vep annotations together if there is more than 
            # one annotation per variant
            vep_annotations = []
            for ann in range(0, len(vep_headers)):
                t = [x[ann] for x in VEP]
                # print(t)
                t = ";".join(t)
                vep_annotations.append(t)
            vep_annotations = vep_annotations + [highest_impact]
            # print("")
        # print(vep_annotations)
        # print("") 
        
        # print(vep_headers)
        # print(sublist_length)
        # print(VEP)
        
        ############################################
        ########## INFO
        
        DP_idx = INFO_names.index('DP')
        INFO_DP = (INFO_values[DP_idx])
        
        MQM_idx = INFO_names.index('MQM')
        INFO_MQM = (INFO_values[MQM_idx])
        
        MQMR_idx = INFO_names.index('MQMR')
        INFO_MQMR = (INFO_values[MQMR_idx])
        
        QA_idx = INFO_names.index('QA')
        INFO_QA = (INFO_values[QA_idx])
        
        QR_idx = INFO_names.index('QR')
        INFO_QR = (INFO_values[QR_idx])
        
        SAF_idx = INFO_names.index('SAF')
        INFO_SAF = (INFO_values[SAF_idx])
        
        SAR_idx = INFO_names.index('SAR')
        INFO_SAR = (INFO_values[SAR_idx])
        
        SRF_idx = INFO_names.index('SRF')
        INFO_SRF = (INFO_values[SRF_idx])
        
        SRR_idx = INFO_names.index('SRR')
        INFO_SRR = (INFO_values[SRR_idx])
        
        SBR_idx = INFO_names.index('SBR')
        INFO_SBR = (INFO_values[SBR_idx])
        
        SBA_idx = INFO_names.index('SBA')
        INFO_SBA = (INFO_values[SBA_idx])
        
        # VAF_idx = INFO_names.index('VAF')
        # INFO_VAF = (INFO_values[VAF_idx])
        # INFO_VAF = INFO_VAF.split(",")
        # print(VAF)
        
        geno_names = lineparts[8]
        geno_names = geno_names.split(":")
        # print(geno_names)
        
        for samp in range(0, n_sample):
            
            # samp_VAF = INFO_VAF[samp]
            # print(samp_VAF)
            
            geno = lineparts[9 + samp]
            geno = geno.split(":")
            # print(geno)
            
            VAF_idx = geno_names.index('VAF')
            FORMAT_VAF = geno[VAF_idx]
            
            if float(FORMAT_VAF) > float(min_vaf):
                
                # print("VAF")
                # print(float(FORMAT_VAF))
                # print(min_vaf)
                # print(float(FORMAT_VAF) > float(min_vaf))
                
                GT_idx = geno_names.index('GT')
                FORMAT_GT = geno[GT_idx]
                # print(GT)
                
                DP_idx = geno_names.index('DP')
                FORMAT_DP = geno[DP_idx]
                # print(FORMAT_DP)
                
                # AD_idx = geno_names.index('AD')
                # FORMAT_AD = str(geno[AD_idx])
                # print(FORMAT_AD)  
                
                RO_idx = geno_names.index('RO')
                FORMAT_RO = geno[RO_idx]
                # print(FORMAT_RO)  
                
                QR_idx = geno_names.index('QR')
                FORMAT_QR = geno[QR_idx]
                # print(FORMAT_QR)  
                
                AQR_idx = geno_names.index('AQR')
                FORMAT_AQR = geno[AQR_idx]
                # print(FORMAT_AQR) 
                
                AO_idx = geno_names.index('AO')
                FORMAT_AO = geno[AO_idx]
                # print(FORMAT_AO)      
                
                QA_idx = geno_names.index('QA')
                FORMAT_QA = geno[QA_idx]
                # print(FORMAT_QA)
                
                AQA_idx = geno_names.index('AQA')
                FORMAT_AQA = geno[AQA_idx]
                # print(FORMAT_AQA)                                          

                # referred to as variant_quality in the final table
                QUAL_idx = geno_names.index('q')
                FORMAT_QUAL = geno[QUAL_idx]
                # print(FORMAT_QUAL)

                # @TODO - move to `mity normalise` as INFO:TIER
                if float(FORMAT_VAF) >= 0.01:
                    tier = 1
                elif float(FORMAT_VAF) < 0.01 and float(FORMAT_AO) > 10:
                    tier = 2
                else:
                    tier = 3
                    
                # table = [["SAMPLE", "CHR", "POS", "REF", "ALT", "QUAL",
                # "FILTER", "INFO_DP", "INFO_MQM", "INFO_MQMR",
                #       "INFO_QA", "INFO_QR", "INFO_SAF", "INFO_SAR",
                #       "INFO_SRF", "INFO_SRR", "INFO_SBR", "INFO_SBA",
                #       "INFO_VAF", "INFO_POS_FILTER", "INFO_SBR_FILTER",
                #       "INFO_SBA_FILTER", "INFO_MQMR_FILTER", "INFO_AQR_FILTER",
                #       "FORMAT_GT", "FORMAT_DP", "FORMAT_RO",
                #       "FORMAT_QR", "FORMAT_AQR", "FORMAT_AO",
                #       "FORMAT_QA", "FORMAT_AQA", "INFO", "FORMAT"]]
                # print(lineparts[7])
                # print(lineparts[9+samp])
                
                no_comma_info = str(lineparts[7])
                no_comma_info = no_comma_info.replace(",", "_")
                
                no_comma_format = str(lineparts[9 + samp])
                no_comma_format = no_comma_format.replace(",", "_")
                
                new_line = [samples[samp], chromosome, pos, ref, alt, hgvs,
                            QUAL, FILTER, INFO_DP, INFO_MQM, INFO_MQMR,
                            INFO_QA, INFO_QR, INFO_SAF, INFO_SAR, INFO_SRF,
                            INFO_SRR, INFO_SBR, INFO_SBA, FORMAT_VAF,
                            POS_FILTER, SBR_FILTER, SBA_FILTER, MQMR_FILTER, AQR_FILTER,
                            FORMAT_GT, FORMAT_DP, FORMAT_RO, FORMAT_QR,
                            FORMAT_AQR, FORMAT_AO,
                            FORMAT_QA, FORMAT_AQA, FORMAT_QUAL, tier, no_comma_info,
                            no_comma_format]
                
                if vep_headers != "":
                    new_line = new_line + vep_annotations

                table.append(new_line)
    
    return (table)

def split_header_variants(vcf):
    header = []
    variants = []
    for line in vcf:
        # print(line)
        if line[0] == '#':
            line = line.strip()
            header.append(line)
            # print(line)
        else:
            line = line.strip()
            # split the line into its columns
            lineparts = line.split('\t')
            variants.append(lineparts)
            
            # print(header)
    # print(variants)
    
    # check if VEP appears anywhere in the header
    is_vepped = ['VEP' in x for x in header]
    is_vepped = True in is_vepped
    # print(is_vepped)
    
    # if its vepped we need to get the vep col names
    if is_vepped:
        csq_line_idx = ['INFO=<ID=CSQ' in x for x in header]
        csq_line_idx = csq_line_idx.index(True)
        csq_line = header[csq_line_idx]
        csq_line = csq_line.strip()
        # remove the last "> from the line)
        vep_header = csq_line[:-2]
        vep_header = vep_header.split("|")
        # remove the first uneeded stuff before Consequence
        vep_header[0] = vep_header[0].rsplit(' ', 1)[1]
    
    
    else:
        vep_header = ""
        # impacts = ""      
        # print(d)
        # print(d['non_coding_transcript_exon_variant'])
    
    # the col names should be the last line of the header
    col_names = header[-1]
    # print(col_names)
    # check that the col names only start with a #
    if not (col_names[0] == '#' and col_names[1] != '#'):
        sys.exit(
            "It looks like the last line of the header is not the column "
            "names.")
    col_names = col_names.strip()
    # remove the first # from CHROM
    col_names = col_names[1:]
    col_names = col_names.split('\t')
    # print(col_names)
    samples = col_names[9:]
    n_sample = len(samples)
    # remove the last n_sample items from the col_names
    del col_names[-n_sample:]
    return (samples, col_names, variants, is_vepped, vep_header, header)


def find_index(string, pattern):
    return [i for i, ltr in enumerate(string) if ltr == pattern]


def do_report(vcf, prefix=None, min_vaf=0.0, out_folder_path = "."):
    """
    Create a mity report
    :param vcf: the path to a vcf file
    :param prefix: the optional prefix. This must be set if there is >1 vcf files
    :param min_vaf: only include vairants with vaf > min_vaf in the report
    :return:
    """
    vcf = vcf[0]
    
    if len(vcf) == 0:
        raise ValueError("At least one VCF file must be supplied")
    if prefix is None and len(vcf) > 1:
        raise ValueError("If there is more than one vcf file, --prefix must be set")
    check_missing_file(vcf, die=True)

    prefix = create_prefix(vcf[0], prefix)
    
    # loop over all the files that are input
    variant_list = []
    
    samples = []
    col_names = []
    variants = []
    is_vepped = []
    vep_header = []
    header_text = []
    
    for filepath in vcf:
        # this assumes they all have the same VEP cols
        file = gzip.open(filepath, 'rt')
        # for each file get the first line that doesnt start with # to get 
        # the sample names
        
        header_variants = split_header_variants(file)
        
        # there will be multiple sample names per vcf (if joint called)
        # each set of samples a sublist. One vcf has a set of samples
        samples.append(header_variants[0])
        
        # one set of col names per vcf
        # each set of column names a sublist
        file_col_names = header_variants[1]
        col_names.append(file_col_names)
        
        # one set of variants per vcf
        # each set of variants a sublist
        file_variants = header_variants[2]
        variants.append(file_variants)
        
        # each vcf is either vepped or not
        # list of true false for "if vepped"
        is_vepped = is_vepped + [header_variants[3]]
        vep_header.append(header_variants[4])
        
        header_text = header_text + header_variants[5]
    # print(header_text)
    # for the Documentation sheet in the output, get the VCF headers
    # we want the header that explains INFO, FILTER or FORMAT fields
    line_types = ["##INFO", "##FILTER", "##FORMAT"]
    vcf_documentation = []
    ID = []
    for line in header_text:
        line_type = line.split("=")[0]
        if line_type in line_types:
            # print(line)
            ##INFO=<ID=DP,Number=1,Type=Integer,Description="
            description_ID = line.split(",")[0]
            ID.append(description_ID)
            description_ID = description_ID.split("ID=")[1]
            # the CSQ header will be taken care of in the VEP columns
            if description_ID != "CSQ":
                line_type = line_type.split("##")[1]
                if line_type != "FILTER":
                    description_ID = description_ID + "_" + line_type + ": "
                else:
                    description_ID = description_ID + ": "
                description = line.split('Description="')[1]
                description = description[:-2]
                new_doc_line = description_ID + description
                vcf_documentation.append(new_doc_line)
    
    noDupes_vcf_documentation = []
    
    for line in vcf_documentation:
        if not noDupes_vcf_documentation.count(line):
            noDupes_vcf_documentation.append(line)
    # print(noDupes_vcf_documentation)
    # if the number of unique IDs does not match the number of unique 
    # descriptions, then there is one ID 
    # that has more than one descrption, either because of one vcf with 
    # duplicated headers, or because
    # 2 VCFs have been inputted with conflicting headers
    dupe_header_warning = False
    if (len(set(ID)) != len(noDupes_vcf_documentation)):
        dupe_header_warning = True
    # print(dupe_header_warning)
    # sys.exit()    
    
    # check they all have the same length
    # print(samples)
    if len(samples) != len(col_names) != len(variants) != len(is_vepped) != len(
            vep_header):
        sys.exit(
            "length of samples, col_names, variants, is_vepped, vep_header is "
            "different, something odd has happened.")
    
    # check that all the col_names are the same
    if len(col_names) > 1:
        l = col_names.pop(0)
        all_equal = [x == l for x in col_names]
        
        if set(all_equal) != set([True]):
            sys.exit(
                "VCF inputs do not have the same column name for the variants")
    
    # if all the col_names are the same we only need to keep one of them
    # print(col_names)
    col_names = col_names[0]
    # print(col_names)
    
    # check they are either all vepped or not vepped
    if len(set(is_vepped)) != 1:
        sys.exit("Some VCFs are vepped and some are not. They must be the same")
    # again, now all we need to know is if they are all vepped or not
    is_vepped = True in is_vepped
    
    # check that the vep headers are the same
    if len(vep_header) > 1:
        l = vep_header.pop(0)
        all_equal = [x == l for x in vep_header]
        
        if set(all_equal) != set([True]):
            sys.exit("VCFs have different VEP columns")
            # if all the vep_header are the same we only need to keep one of 
            # them
    vep_header = vep_header[0]
    
    # if its vepped we need the impacts
    if is_vepped:
        impacts = {}
        # http://asia.ensembl.org/info/genome/variation/predicted_data.html
        with open("ensemble_impacts.txt") as f:
            for line in f:
                (key, val) = line.split()
                impacts[key] = val
    else:
        impacts = ""
    
    # range!!!!!
    
    ## get the variants into a list format
    variant_table = []
    # print(len(variants))
    # print("making table")
    for vcf_num in range(0, len(variants)):
        v = make_table(variants[vcf_num], samples[vcf_num], vep_header, impacts,
                       min_vaf)
        # variant_table.append(v)w
        variant_table = variant_table + v
    # To print just the variants:
    # for line in variant_table:
    # print(line)
    # sys.exit()    
    
    # note that this is hard coded and the function make_table outputs the 
    # columns in this order
    # the column order is changed later in the code
    vcf_headers = ["SAMPLE", "CHR", "POS", "REF", "ALT", "HGVS", "QUAL",
                   "FILTER", "total_locus_depth", "MQM_INFO", "MQMR_INFO",
                   "QA_INFO", "QR_INFO", "SAF_INFO", "SAR_INFO", "SRF_INFO",
                   "SRR_INFO", "SBR_INFO", "SBA_INFO",
                   "variant_heteroplasmy", "POS_FILTER", "SBR_FILTER", "SBA_FILTER",
                   "MQMR_FILTER", "AQR_FILTER",
                   "GT_FORMAT", "total_sample_depth", "ref_depth", "QR_FORMAT",
                   "AQR_FORMAT", "alt_depth",
                   "QA_FORMAT", "AQA_FORMAT", "variant_quality", "tier", "INFO", "FORMAT"]
    
    if vep_header != "":
        table_vep_headers = [x + "_VEP" for x in vep_header]
        # we add 'highest_vep_impact' to the end of the VEP cols
        vcf_headers = vcf_headers + table_vep_headers
        vcf_headers = vcf_headers + ['highest_vep_impact']
        # print(headers)
    # print(headers)
    
    variant_df = pandas.DataFrame(variant_table, columns=vcf_headers)
    # print(variant_df)
    # variant_df.to_csv("test_variants.csv", index = False)

    # The annotations use contig name = 'MT', so if the vcf uses contig name = 'chrM', then swap this.
    contig_name = variant_df['CHR'][0]
    if contig_name != "MT":
        variant_df['CHR'] = "MT"
    ########## Merge with mitomap and panel annotation table
    # this depends on chrom, pos, ref, alt
    mitomap_panel_annotations = pandas.read_csv(get_annot_file("mitomap_panel_annotations.csv"))
    mitomap_panel_annotations['POS'] = mitomap_panel_annotations['POS'].astype('str')
    # print(mitomap_panel_annotations.dtypes)
    # print(variant_df.dtypes)
    mitomap_panel_annotated_variants = pandas.merge(left=variant_df,
                                                    right=mitomap_panel_annotations,
                                                    how='left',
                                                    on=['CHR', 'POS', 'REF', 'ALT'])
    
    ########## Merge with gene names and biotype
    # this depends on chrom, pos
    gtf_annotations = pandas.read_csv(get_annot_file("gtf_annotations.csv"))
    gtf_annotations['POS'] = gtf_annotations['POS'].astype('str')
    # print(gtf_annotations.dtypes)
    
    gtf_mitomap_panel_annotated_variants = pandas.merge(
        left=mitomap_panel_annotated_variants, right=gtf_annotations,
        how='left', on=['CHR', 'POS'])
    # annotated_variants = annotated_variants.fillna('.')
    
    ########## Merge with Mitomap Locations
    # this depends on chrom, pos
    mito_locus_annotations = pandas.read_csv(get_annot_file("mito_dna_func_loc.csv"))
    mito_locus_annotations['POS'] = mito_locus_annotations['POS'].astype('str')
    # print(gtf_annotations.dtypes)
    
    mito_locus_gtf_mitomap_panel_annotated_variants = pandas.merge(
        left=gtf_mitomap_panel_annotated_variants, right=mito_locus_annotations,
        how='left', on=['CHR', 'POS'])
    # for line in mito_locus_gtf_mitomap_panel_annotated_variants:
    #   print(line)
    # sys.exit()
    
    ########## Merge with trna anticodon positions
    # this depends on chrom pos 
    trna_annotations = pandas.read_csv(get_annot_file("anticodon_positions.csv"))
    trna_annotations['POS'] = trna_annotations['POS'].astype('str')
    trna_annotations['anticodon'] = trna_annotations['anticodon'].astype('str')
    trna_mito_locus_gtf_mitomap_panel_annotated_variants = pandas.merge(
        left=mito_locus_gtf_mitomap_panel_annotated_variants,
        right=trna_annotations, how='left', on=['CHR', 'POS'])
    # for line in trna_mito_locus_gtf_mitomap_panel_annotated_variants:
    #   print(line)
    # sys.exit()
    
    ########## Merge with mitotip score and prediction
    # this depends on chrom pos ref alt
    mitotip_annotations = pandas.read_csv(get_annot_file("mitotip_score_fixed_del.csv"))
    mitotip_annotations['POS'] = mitotip_annotations['POS'].astype('str')
    mitotip_annotations['MitoTip_score'] = mitotip_annotations[
        'MitoTip_score'].astype('str')
    mitotip_annotations['MitoTip_percentile'] = mitotip_annotations[
        'MitoTip_percentile'].astype('str')
    mitotip_annotations['MitoTip_interpretation'] = mitotip_annotations[
        'MitoTip_interpretation'].astype('str')
    mitotip_trna_mito_locus_gtf_mitomap_panel_annotated_variants = pandas.merge(
        left=trna_mito_locus_gtf_mitomap_panel_annotated_variants,
        right=mitotip_annotations, how='left', on=['CHR', 'POS', 'REF', 'ALT'])
    # for line in trna_mito_locus_gtf_mitomap_panel_annotated_variants:
    #   print(line)
    # sys.exit()
    
    ########## Merge with mgrb variants
    # this depends on chrom pos ref alt
    mgrb_annotations = pandas.read_csv(get_annot_file("mgrb_variants.csv"))
    mgrb_annotations['POS'] = mgrb_annotations['POS'].astype('str')
    mgrb_mitotip_trna_mito_locus_gtf_mitomap_panel_annotated_variants = \
        pandas.merge(
        left=mitotip_trna_mito_locus_gtf_mitomap_panel_annotated_variants,
        right=mgrb_annotations, how='left', on=['CHR', 'POS', 'REF', 'ALT'])
    # for line in trna_mito_locus_gtf_mitomap_panel_annotated_variants:
    # print(line)
    # sys.exit()
    
    ########## Merge with haplotype data
    # TODO: why does this add two lines for some variants? (there are no pos-ref-alt dups)
    haplotype_annotations = pandas.read_csv(get_annot_file("haplotype_data.csv"))
    haplotype_annotations['POS'] = haplotype_annotations['POS'].astype('str')
    haplotype_mgrb_mitotip_trna_mito_locus_gtf_mitomap_panel_annotated_variants = pandas.merge(
        left=mgrb_mitotip_trna_mito_locus_gtf_mitomap_panel_annotated_variants,
        right=haplotype_annotations, how='left',
        on=['CHR', 'POS', 'REF', 'ALT'])
    # for line in trna_mito_locus_gtf_mitomap_panel_annotated_variants:
    # print(line)
    # sys.exit()
    
    annotated_variants = \
        haplotype_mgrb_mitotip_trna_mito_locus_gtf_mitomap_panel_annotated_variants.fillna(
        '.')
    # print(annotated_variants)
    # sys.exit()

    # The annotations use contig name = 'MT', so if the vcf uses contig name = 'chrM', then swap this.
    if contig_name != "MT":
        annotated_variants["CHR"] = contig_name

    # now we make the general gene/locus and gene/locus description cols
    gene = annotated_variants['GENE'].tolist()
    gene_biotype = annotated_variants['GENE_BIOTYPE'].tolist()
    mitomap_locus = annotated_variants['Map_Locus'].tolist()
    mitomap_description = annotated_variants['Description'].tolist()
    
    general_locus = gene
    general_biotype = gene_biotype
    # now replace "."'s with mitomap entries
    empty_idx = find_index(general_locus, ".")
    for x in empty_idx:
        general_locus[x] = mitomap_locus[x]
        general_biotype[x] = mitomap_description[x]
    
    # we remove gene, gene biotype, Map_Locus and Description from table
    del annotated_variants['GENE']
    del annotated_variants['GENE_BIOTYPE']
    del annotated_variants['Map_Locus']
    del annotated_variants['Description']
    del annotated_variants['baylor_panel']
    del annotated_variants['common_22_panel']
    del annotated_variants['common_58_panel']
    #annotated_variants1.drop(columns=['baylor_panel','common_22_panel','common_58_panel'])

    # add the general locus and general locus description
    annotated_variants['gene/locus'] = general_locus
    annotated_variants['gene/locus description'] = general_biotype
    
    cohort_count_table = pandas.DataFrame(
        {'count': annotated_variants.groupby(["POS", "ALT"]).size()}).reset_index()
    cohort_count_table.columns = ['POS', 'ALT', 'COHORT_COUNT']
    
    annotated_variants1 = pandas.merge(left=annotated_variants,
                                       right=cohort_count_table, how='left',
                                       on=['POS', 'ALT'])
    
    # add in the GB percent
    # need to account for the fact that some entries will be "." and cant be 
    # converted to float
    gb_freq = annotated_variants1['GenBank_frequency_mitomap'].tolist()
    gb_perc = gb_freq
    
    # only convert the entries that are not a "."
    str_idx = find_index(gb_freq, ".")
    str_idx = [i for i, ltr in enumerate(gb_freq) if ltr != "."]
    
    for x in str_idx:
        freq = round(float(gb_freq[x]) / 32059, 3)
        # print(freq)
        gb_perc[x] = freq
    
    annotated_variants1['allele_frequency_mitomap'] = gb_perc
    # TODO: to make this more general, choose some important cols to go first
    # if they are in the header then do this sorting
    # otherwise dont
    if is_vepped:
        cols = ['SAMPLE', 'HGVS', 'gene/locus', 'gene/locus description',
                'variant_heteroplasmy', 'alt_depth',
                'ref_depth', 'total_sample_depth', 'variant_quality',
                'total_locus_depth', 'COHORT_COUNT', 'tier', 'commercial_panels',
                'phylotree_haplotype', 'MitoTip_score',
                'MitoTip_percentile', 'MitoTip_interpretation', 'anticodon',
                'allele_frequency_mitomap',
                'highest_vep_impact', 'Consequence_VEP', 'disease_mitomap',
                'MGRB_frequency', 'MGRB_FILTER', 'MGRB_AC', 'MGRB_AN', 'phylotree_mut',
                'Codons_VEP', 'Amino_acids_VEP', 'Gene_VEP', 'SYMBOL_VEP',
                'Feature_VEP', 'EXON_VEP', 'PolyPhen_VEP', 'SIFT_VEP',
                'Protein_position_VEP', 'BIOTYPE_VEP', 'CANONICAL_VEP',
                'Feature_type_VEP', 'cDNA_position_VEP', 'CDS_position_VEP',
                'Existing_variation_VEP', 'DISTANCE_VEP', 'STRAND_VEP',
                'CLIN_SIG_VEP', 'LoF_flags_VEP', 'LoF_filter_VEP', 'LoF_VEP',
                'RadialSVM_score_VEP', 'RadialSVM_pred_VEP', 'LR_score_VEP',
                'LR_pred_VEP', 'CADD_raw_VEP', 'CADD_phred_VEP',
                'Reliability_index_VEP', 'HGVSc_VEP', 'HGVSp_VEP',
                'locus_mitomap', 'num_references_mitomap',
                'variant_amino_acid_change_mitomap', 'codon_position_mitomap',
                'codon_number_mitomap',
                'num_disease_references_mitomap', 'RNA_mitomap',
                'homoplasmy_mitomap', 'heteroplasmy_mitomap', 'status_mitomap',
                'disease_amino_acid_change_mitomap', 'GenBank_frequency_mitomap',
                'CHR', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER',
                'MQM_INFO', 'MQMR_INFO', 'QA_INFO', 'QR_INFO', 'SAF_INFO',
                'SAR_INFO', 'SRF_INFO', 'SRR_INFO', 'SBR_INFO', 'SBA_INFO',
                'POS_FILTER', 'SBR_FILTER', 'SBA_FILTER', 'MQMR_FILTER', 'AQR_FILTER',
                'GT_FORMAT', 'QR_FORMAT', 'AQR_FORMAT', 'QA_FORMAT',
                'AQA_FORMAT',
                'INFO', 'FORMAT']
    else:
        cols = ['SAMPLE', 'HGVS', 'gene/locus', 'gene/locus description',
                'variant_heteroplasmy', 'alt_depth',
                'ref_depth', 'total_sample_depth', 'variant_quality',
                'total_locus_depth', 'COHORT_COUNT', 'tier', 'commercial_panels',
                'phylotree_haplotype', 'MitoTip_score',
                'MitoTip_percentile', 'MitoTip_interpretation', 'anticodon',
                'allele_frequency_mitomap',
                'disease_mitomap',
                'MGRB_frequency', 'MGRB_FILTER', 'MGRB_AC', 'MGRB_AN',
                'phylotree_mut',
                'locus_mitomap', 'num_references_mitomap',
                'variant_amino_acid_change_mitomap', 'codon_position_mitomap',
                'codon_number_mitomap',
                'num_disease_references_mitomap', 'RNA_mitomap',
                'homoplasmy_mitomap', 'heteroplasmy_mitomap', 'status_mitomap',
                'disease_amino_acid_change_mitomap', 'GenBank_frequency_mitomap',
                'CHR', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER',
                'MQM_INFO', 'MQMR_INFO', 'QA_INFO', 'QR_INFO', 'SAF_INFO',
                'SAR_INFO', 'SRF_INFO', 'SRR_INFO', 'SBR_INFO', 'SBA_INFO',
                'POS_FILTER', 'SBR_FILTER', 'SBA_FILTER', 'MQMR_FILTER', 'AQR_FILTER',
                'GT_FORMAT', 'QR_FORMAT', 'AQR_FORMAT', 'QA_FORMAT',
                'AQA_FORMAT',
                'INFO', 'FORMAT']
        # These should be equal if you haven't missed any cols
    
    if len(annotated_variants1.columns.tolist()) != len(cols):
        print(cols)
        print("")
        print(annotated_variants.columns.tolist())
        print("")
        print("annotated_variants missing" + str(
            set(annotated_variants.columns.tolist()) - set(cols)))
        print("")
        print(
            "cols" + str(set(cols) - set(annotated_variants.columns.tolist())))
        sys.exit("different number of cols")
    
    annotated_variants1 = annotated_variants1[cols]
    
    # make the documentation sheet
    documentation = []
    x = "This sheet contains the variants from"
    for filename in vcf:
        name = filename.split("/")[-1]
        x = x + " " + str(name)
    documentation.append(x)
    documentation.append("Some column definitions:")
    documentation.append(
        "HGVS: Variant annotation using HGVS syntax, which assumes variant is "
        "SNP, insertion or deletion")
    documentation.append(
        "gene/locus: The gene that the variant falls in, or if it is not in a "
        "gene then the mitomap entry, if it exists.")
    documentation.append(
        "variant_heteroplasmy: alt_depth/(ref_depth+alt_depth)")
    documentation.append(
        "variant_quality: phred-scaled variant quality")
    documentation.append(
        "total_sample_depth: alt_depth + ref_depth for the particular sample")
    documentation.append(
        "total_locus_depth: sum of all the total_sample_depth in the vcf file")
    documentation.append(
        "COHORT_COUNT: The number of times that variant (defined by position "
        "and alternate) is seen in all of the input VCFs")
    documentation.append(
        "tier: 1:=VAF>=1%. 2:= VAF<1% and alt_depth>=10, 3:=VAF<1% and "
        "alt_depth<10")
    documentation.append(
        "commercial_panels: This indicates whether a variant is commonly tested "
        "in a number of commercial/academic testing laboratories. "
        "Baylor: https://www.bcm.edu/research/medical-genetics-labs/test_detail.cfm?testcode=2010; "
        "Common22: https://www.vcgs.org.au/order/tests/634; "
        "Common58: https://www.genedx.com/test-catalog/available-tests"
        "/58-confirmed-disease-causing-mtdna-point-mutations-and-deletion"
        "-testing/")
    documentation.append(
        "MGRB: The Medical Genomics Reference Bank (MGRB) "
         "is a cohort of healthy older Ausrtralian's sequenced using "
         "~30x Illumina HiSeq X at the Garvan Institute, Sydney. "
         "See https://www.biorxiv.org/content/10.1101/473348v1. ")
    documentation.append(
        "MGRB_frequency: represents the frequency of this variant in the "
        "MGRB cohort")
    documentation.append(
        "MGRB_FILTER: represents the FILTER status of the variant in the "
        "MGRB cohort using GATK HaplotypeCaller")
    documentation.append(
        "MGRB_AC: the allele count in the MGRB cohort; ie the number of "
        "individuals with this variant. "
        "This is not accurate for the MT, as it assumes each individual is "
        "diploid. Interpret this value as 2x the number of individuals "
        "in the MGRB control cohort.")
    documentation.append(
        "MGRB_AN: the total number of alleles called in the MGRB cohort. "
        "This is not accurate for the MT, as it assumes each individual is "
        "diploid. Interpret this value as 2x the number of individuals "
        "in the MGRB control cohort.")
    documentation.append(
        "phylotree_haplotype: The haplotypes that the variant is known "
        "to contribute to. See http://www.phylotree.org/index.htm")
    documentation.append(
        "MitoTip_score: MitoTip is an in silico tool for predicting the "
        "pathogenicity of novel mitochondrial tRNA variants, developed "
        "by Neal Sondheimer and Sanjay Sonney (PMID: 29227991). "
        "As such, it only scores variants in a tRNA. The score is a "
        "prediction of the pathogenicity of variants in tRNAs, "
        "see https://www.mitomap.org/foswiki/bin/view/MITOMAP/MitoTipInfo. "
        "The highest score is 21.8. It is far from perfect, as m.3243A>G "
        "is scored as likely benign. ")
    documentation.append(
        "MitoTip_percentile: Where this variants MitoTip score sits in "
        "relation to all of the possible variants in tRNA. 100% is the most "
        "pathogenic, 0% is a score of zero.")
    documentation.append(
        "MitoTip_interpretation: Variants in the upper 25th percentile are "
        "deemed likely pathogenic, 50-75th percentile as possibly "
        "pathogenic; 25-50th percentile as possibly benign; "
        "bottom 25th percentile as likely benign. "
        "Variants reported as Confirmed pathogenic in MITOMAP have additional lines " 
        "of experimental evidence which led to confirmation of pathogenicity "
        "regardless of which percentile category the variant may belong to. "
        "see "
        "https://mitomap.org/foswiki/bin/view/MITOMAP/MitoTipInfo")
    documentation.append(
        "anticodon: TRUE if the variant falls in an anticodon of a tRNA.")
    documentation.append(
        "Columns suffix with _mitomap are annotations drawn from the MITOMAP "
        "compendium of polymorphisms and mutations in human mitochondrial DNA. "
        "see http://www.mitomap.org/MITOMAP. "
        "Variants from the the 4 tables are included: 'Control Region Variants "
        "(16024-576)', 'Coding & RNA Variants (577-16023, MTTF-MTTP)', 'rRNA/tRNA "
        " Mutations' and 'Coding & Control Region Mutations' (obtained in 2018)")
    documentation.append(
        "disease_mitomap: The disease annotation for a particular variant, "
        "as annotated by MITOMAP")
    documentation.append(
        "status_mitomap: an indication of the strength of evidence supporting "
        "the MITOMAP annotation. The strongest evidence is Confirmed, which "
        "MITMAP calls Cfrm. "
        "See https://mitomap.org/foswiki/bin/view/MITOMAP/MutationsRNACfrm")
    documentation.append(
        "GenBank_frequency_mitomap: the frequency of GenBank records "
        "containing this variant. This is a proxy for the population "
        "frequency of the variant. "
        "As expected, none of the confirmed pathogenic variants in MITOMAP "
        "have GenBank frequency >0.0%")
    documentation.append(
        "allele_frequency_mitomap: MITOMAP derives this from 32,059 GenBank "
        "full length sequences, with size greater than 15.4kbp (range 0-1)")
    documentation.append(
        "Columns suffix with _VEP come from Ensembl’s Variant Effect "
        "Predictor (VEP) annotations in the VCF")
    documentation.append(
        "highest_vep_impact: Highest VEP impact from all the VEP annotations "
        "for the variant. See http://www.ensembl.org/info/docs/tools/vep/index.html")
    documentation.append(
        "Entries separated with a semicolon mean that there are multiple "
        "annotations, from either VEP or MITOMAP.")
    documentation.append(
        "Fields from VCF header: CHR, POS, REF, ALT, QUAL, FILTER, "
        "INFO, FORMAT")
    documentation.append(
        "Additional fields from mity VCF INFO field: MQM_INFO, "
        "MQMR_INFO, QA_INFO, QR_INFO, SAF_INFO, SAR_INFO, SRF_INFO, "
        "SRR_INFO, SBR_INFO, SBA_INFO")
    documentation.append(
        "Additional fields from the mity VCF FILTER field: POS_FILTER, "
        "SBR_FILTER, SBA_FILTER, MQMR_FILTER, AQR_FILTER")
    documentation.append(
        "Additional fields from the mity VCF FORMAT field: GT_FORMAT, "
        "QR_FORMAT, AQR_FORMAT, QA_FORMAT, AQA_FORMAT")
    if dupe_header_warning:
        documentation.append(
            "Warning: VCF header IDs are not unique, so there might be a "
            "repeated definition below.")
    
    for line in noDupes_vcf_documentation:
        documentation.append(line)
    
    documentation_df = pandas.DataFrame(documentation)
    
    ## change the columns that are numerical data types to numerical
    annotated_variants1['variant_heteroplasmy'] = annotated_variants1[
        'variant_heteroplasmy'].astype('float64')
    annotated_variants1['variant_quality'] = annotated_variants1[
        'variant_quality'].astype('float64')
    annotated_variants1['ref_depth'] = annotated_variants1['ref_depth'].astype(
        'int64')
    annotated_variants1['alt_depth'] = annotated_variants1['alt_depth'].astype(
        'int64')
    annotated_variants1['total_sample_depth'] = annotated_variants1[
        'total_sample_depth'].astype('int64')
    annotated_variants1['total_locus_depth'] = annotated_variants1[
        'total_locus_depth'].astype('int64')
    annotated_variants1['COHORT_COUNT'] = annotated_variants1[
        'COHORT_COUNT'].astype('int64')
    annotated_variants1['tier'] = annotated_variants1['tier'].astype('int64')
    annotated_variants1['MQM_INFO'] = annotated_variants1['MQM_INFO'].astype(
        'float64')
    annotated_variants1['MQMR_INFO'] = annotated_variants1['MQMR_INFO'].astype(
        'float64')
    annotated_variants1['QA_INFO'] = annotated_variants1['QA_INFO'].astype(
        'float64')
    annotated_variants1['QR_INFO'] = annotated_variants1['QR_INFO'].astype(
        'float64')
    annotated_variants1['SAF_INFO'] = annotated_variants1['SAF_INFO'].astype(
        'float64')
    annotated_variants1['SAR_INFO'] = annotated_variants1['SAR_INFO'].astype(
        'float64')
    annotated_variants1['SRF_INFO'] = annotated_variants1['SRF_INFO'].astype(
        'float64')
    annotated_variants1['SRR_INFO'] = annotated_variants1['SRR_INFO'].astype(
        'float64')
    annotated_variants1['SBR_INFO'] = annotated_variants1['SBR_INFO'].astype(
        'float64')
    annotated_variants1['SBA_INFO'] = annotated_variants1['SBA_INFO'].astype(
        'float64')

    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    xlsx_name = os.path.join(out_folder_path, prefix + ".annotated_variants.xlsx")
    logging.info("saving xlsx report: " + xlsx_name)
    writer = pandas.ExcelWriter(xlsx_name, engine='xlsxwriter')
    documentation_df.to_excel(writer, sheet_name='Documentation', index=False,
                              header=False)
    annotated_variants1.to_excel(writer, sheet_name='Variants', index=False)
    writer.save()
    
    csv_name = os.path.join(out_folder_path, prefix + ".annotated_variants.csv")
    logging.info("saving csv report: " + csv_name)
    annotated_variants1.to_csv(csv_name, index=False, )
    # numpy.savetxt('data.csv', delimiter=',', X = annotated_variants1)
    # with open("test.csv", "w", newline='') as csv_file:
    #   writer = csv.writer(csv_file, delimiter=',')
    #   for line in annotated_variants1:
    #       writer.writerow(line)
