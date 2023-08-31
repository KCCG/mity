"""Mitochondrial variant calling."""
import sys
import subprocess
import logging
import os.path
from .util import tabix, check_missing_file, create_prefix, bam_get_mt_contig, bam_has_RG
from .normalise import do_normalise as vcfnorm

# CONSTANTS
MIN_MQ=30
MIN_BQ=24
MIN_AF=0.01
MIN_AC=4
P_VAL=0.002

# LOGGING
logger = logging.getLogger(__name__)

def do_call(debug, bam_files, reference, genome=None, prefix=None, min_mq=MIN_MQ, min_bq=MIN_BQ,
            min_af=MIN_AF, min_ac=MIN_AC, p=P_VAL, normalise=True, 
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

    if debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Entered debug mode.")
    else:
        logger.setLevel(logging.INFO)

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
        logger.error("At least one BAM/CRAM file lacks an @RG header")
        exit(1)

    if normalise and genome is None:
        logger.error("A genome file should be supplied if mity call normalize=True")
        exit(1)

    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    mity_normalise_output_file = os.path.join(out_folder_path, prefix + ".mity.normalise.vcf.gz")
    mity_call_output_file = os.path.join(out_folder_path, prefix + ".mity.call.vcf.gz")

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
    logger.debug("mity commandline: " + str(mity_cmd))

    if normalise:
        mity_cmd = mity_cmd + ' --normalise --p ' + str(p)

    mity_cmd = mity_cmd + ' ' + ' '.join(bam_files)

    mity_cmd = mity_cmd + '"'
    mity_cmd = mity_cmd.replace("/", "\/")
    logger.debug(mity_cmd)

    # overwrite a redundant freebayes header line with the mity command line
    sed_cmd = "sed 's/^##phasing=none/{}/g'".format(mity_cmd)
    logger.debug(sed_cmd)

    freebayes_call = ('set -o pipefail && freebayes -f {} {} '
                      '--min-mapping-quality {} '
                      '--min-base-quality {} '
                      '--min-alternate-fraction {} '
                      '--min-alternate-count {} '
                      '--ploidy 2 '
                      '--region {} '
                      ).format(reference, bam_str, min_mq, min_bq, min_af, min_ac, region)
    freebayes_call = freebayes_call + ('| sed "s/##source/##freebayesSource/" | sed "s/##commandline/##freebayesCommandline/" | {} | bgzip > {} ').format(sed_cmd, mity_call_output_file)

    logger.info("Running FreeBayes in sensitive mode")
    logger.debug(freebayes_call)
    res = subprocess.run(freebayes_call, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, executable='/bin/bash')
    logger.debug("Freebayes result code: {}".format(res.returncode))

    if res.returncode != 0:
        logger.error("FreeBayes failed: {}".format(res.stderr))
        exit(1)

    if os.path.isfile(mity_call_output_file):
        logger.debug("Finished running FreeBayes")
    
    if normalise:
        logger.debug("Normalising and Filtering variants")

        # Sets the normalise debug flag to true if mity call was called with --debug
        debug_normalise = False
        if debug:
            debug_normalise = True

        try:
            vcfnorm(debug_normalise, vcf=mity_call_output_file, reference_fasta=reference, output_file=mity_normalise_output_file, allsamples=False, p=p, genome=genome)
        finally:
            os.remove(mity_call_output_file)
    else:
        tabix(mity_call_output_file)
