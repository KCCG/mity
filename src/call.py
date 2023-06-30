"""Mitochondrial variant calling."""
import sys
import subprocess
import logging
import os.path
from .util import tabix, check_missing_file, create_prefix, bam_get_mt_contig, bam_has_RG
from .normalise import do_normalise as vcfnorm

def do_call(bam_files, reference, genome=None, prefix=None, min_mq=30, min_bq=24,
            min_af=0.01, min_ac=4, p=0.002, normalise=True, 
            out_folder_path=".", region=None):
    """
    Run mity call.
    :param bam_files: a list of bam or cram files
    :param reference: the path to the reference genome file (fasta format)
    :param genome: the path to the reference genome file for gsort (genome format). Required if normalise=True
    :param prefix: The result filename prefix. If None, then the first bam_file prefix
    will be used
    :param min_mq: minimum mapping quality threshold. default: 30
    :param min_bq: minimum base quality threshold. default: 24
    :param min_af: minimum heteroplasmy, aka the minimum fraction of alt reads vs total reads.
    scale [0,1]; default 0.005
    :param min_ac: minimum number of alternative reads to support a variant. default: 4
    :param p: the noise threshold. default 0.002
    :param normalise:
    :param out_folder_path: the folder to store the results within. default: .
    :param region: Which region to analyse? If None, then the whole MT will be analysed.
    :return:
    """
    bam_files = bam_files[0] 
    #####
    # Checks
    #####

    if len(bam_files) > 1 and prefix is None:
        raise ValueError(
                "If there is more than one bam/cram file, --prefix must be set")
    
    check_missing_file(bam_files, die=True)
    prefix = create_prefix(bam_files[0], prefix)

    if not all(map(bam_has_RG, bam_files)):
        logging.error("At least one BAM/CRAM file lacks an @RG header")
        exit(1)

    if normalise and genome is None:
        logging.error("A genome file should be supplied if mity call normalize=True")
        exit(1)

    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    output_file_name = os.path.join(out_folder_path, prefix + ".mity.vcf.gz")
    unnormalised_vcf_path = os.path.join(out_folder_path, prefix + ".unnormalised.vcf.gz")

    # if given BAM files A, B, C; freebayes will put samples C, B, A in the VCF file; reverse this odd behaviour
    bam_str = " ".join(['-b ' + bam_file for bam_file in reversed(bam_files)])

    if region is None:
        region = bam_get_mt_contig(bam_files[0], as_string=True)

    # embed the mity command into the VCF header
    mity_cmd = '##mityCommandline="mity call --reference ' + str(reference) + \
               ' --prefix ' + prefix + ' --min-mapping-quality ' + str(min_mq) + \
               ' --min-base-quality ' + str(min_bq) + ' --min-alternate-fraction ' + \
               str(min_af) + ' --min-alternate-count ' + str(min_ac) + \
               ' --out-folder-path ' + str(out_folder_path) + ' --region ' + region
    logging.debug("mity commandline: " + str(mity_cmd))

    if normalise:
        mity_cmd = mity_cmd + ' --normalise --p ' + str(p)

    mity_cmd = mity_cmd + ' ' + ' '.join(bam_files)

    mity_cmd = mity_cmd + '"'
    mity_cmd = mity_cmd.replace("/", "\/")
    logging.debug(mity_cmd)

    # overwrite a redundant freebayes header line with the mity command line
    sed_cmd = "sed 's/^##phasing=none/{}/g'".format(mity_cmd)
    logging.debug(sed_cmd)

    freebayes_call = ('set -o pipefail && freebayes -f {} {} '
                      '--min-mapping-quality {} '
                      '--min-base-quality {} '
                      '--min-alternate-fraction {} '
                      '--min-alternate-count {} '
                      '--ploidy 2 '
                      '--region {} '
                      ).format(reference, bam_str, min_mq, min_bq, min_af, min_ac, region)
    freebayes_call = freebayes_call + ('| sed "s/##source/##freebayesSource/" | sed "s/##commandline/##freebayesCommandline/" | {} | bgzip > {} ').format(sed_cmd, unnormalised_vcf_path)

    logging.info("Running FreeBayes in sensitive mode")
    logging.debug(freebayes_call)
    res = subprocess.run(freebayes_call, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable='/bin/bash')
    logging.debug("Freebayes result code: {}".format(res.returncode))

    if res.returncode != 0:
        logging.error("FreeBayes failed: {}".format(res.stderr))
        exit(1)

    if os.path.isfile(unnormalised_vcf_path):
        logging.debug("Finished running FreeBayes")
    
    if normalise:
        logging.debug("Normalising and Filtering variants")
        try:
            vcfnorm(vcf=unnormalised_vcf_path, out_file=output_file_name, p=p, genome=genome)
        finally:
            os.remove(unnormalised_vcf_path)
    else:
        os.rename(unnormalised_vcf_path, output_file_name)
        tabix(output_file_name)
