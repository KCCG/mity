import configparser
import logging
import os
import subprocess
import sys
import tempfile
import vcf as pyvcf
import pysam
import inspect
from glob import glob

def tabix(f):
    """
    Generate a tabix index for a bgzipped file

    :param f: path to a bgzip compressed file
    :type f: str

    :returns: Nothing
    :rtype: None
    """

    tabix_call = "tabix -f " + f
    logging.debug(tabix_call)
    subprocess.run(tabix_call, shell=True)

def check_missing_file(file_list, die=True):
    missing_files = []
    for item in file_list:
        if not os.path.isfile(item):
            missing_files.append(item)
    if die and len(missing_files) > 0:
        raise ValueError("Missing these files: " + ",".join(missing_files))
    return missing_files

def tmp_mity_file_name():
    """
    Create a tmp mity vcf file

    #TODO There must be a more pythonic way of doing this.
    """
    f = tempfile.NamedTemporaryFile(mode="wt", prefix='mity', suffix=".vcf",
                                    delete=False)
    f.close()
    assert isinstance(f.name, str)
    return f.name


def create_prefix(file_name, prefix=None):
    """
    Most mity functions have an optional prefix. If a prefix is not specified,
    then use the  file name (minus the .vcf.gz, .bam or .cram suffix) as the prefix.
    :param file_name: The vcf, bam, cram, bed filename
    :param prefix: The optional prefix. If None, then create a prefix from
    file_name, else return prefix
    :return: str prefix
    """
    if prefix is not None:
        pass
    elif ".vcf" in file_name:
        prefix = [os.path.basename(file_name).split(".vcf")[0], prefix][prefix is not None]
    elif ".bam" in file_name:
        prefix = [os.path.basename(file_name).split(".bam")[0], prefix][prefix is not None]
    elif ".cram" in file_name:
        prefix = [os.path.basename(file_name).split(".cram")[0], prefix][prefix is not None]
    else:
        raise ValueError("Unsupported file type")
    return prefix


def write_vcf(new_vcf, out_file, genome_file='mitylib/reference/b37d5.genome'):
    """
    write a vcf object to vcf.gz file with tbi index.

    This differs from write_merged_vcf, as mity merge doesn't split each VCF
    line on tabs, whereas mity normalise does

    :param new_vcf: new_vcf is a list of lists, created by normalise
    :param out_file: the resulting filename. this should end in vcf.gz
    :return: None. This function writes a vcf.gz and vcf.gz.tbi file.
    """
    f = tmp_mity_file_name()
    logging.debug("Writing uncompressed vcf to " + f)
    with open(f, mode='wt', encoding='utf-8') as myfile:
        for vcf_line in new_vcf:
            myfile.write('\t'.join([str(elem) for elem in vcf_line]) + '\n')
    gsort_vcf(f, out_file, genome_file=genome_file)
    # bcftools_sort_vcf(f, out_file)


def gsort_vcf(f, out_file, genome_file='mitylib/reference/b37d5.genome', remove_unsorted_vcf=False):
    """
    use gsort to sort the records in a VCF file according to a .genome file.

    :param f: the path to an unsorted vcf.gz file
    :param out_file: the path to a resulting sorted vcf.gz file
    :param genome_file: the .genome file corresponding to the reference genome. see https://github.com/brentp/gsort
    :param remove_unsorted_vcf: if True, then the input file 'f' will be deleted.
    :return: nothing
    """
    logging.debug("Sorting, bgzipping {} -> {}".format(f, out_file))
    logging.debug("gsort is using genome file " + genome_file)
    gsort_cmd = "gsort {} {} | bgzip -cf > {}".format(f, genome_file, out_file)

    logging.debug(gsort_cmd)

    subprocess.run(gsort_cmd, shell=True)
    logging.debug("Tabix indexing {}".format(out_file))
    tabix(out_file)
    if remove_unsorted_vcf:
        os.remove(f)

def bcftools_sort_vcf(f, out_file, remove_unsorted_vcf=False):
    """
    use bcftools sort to sort the records in a VCF file.

    :param f: the path to an unsorted vcf.gz file
    :param out_file: the path to a resulting sorted vcf.gz file
    :param remove_unsorted_vcf: if True, then the input file 'f' will be deleted.
    :return: nothing
    """
    logging.debug("Sorting, bgzipping {} -> {}".format(f, out_file))
    bcftools_sort_cmd = "bcftools sort {} -O z -o {} 2>&1 >/dev/null".format(f, out_file)

    logging.debug(bcftools_sort_cmd)

    subprocess.run(bcftools_sort_cmd, shell=True)
    logging.debug("Tabix indexing {}".format(out_file))
    tabix(out_file)
    if remove_unsorted_vcf:
        os.remove(f)

def write_merged_vcf(new_vcf, out_file, genome_file='mitylib/reference/b37d5.genome'):
    """
    write a vcf object to vcf.gz file with tbi index.

    This differs from write_vcf, as mity merge doesn't split each VCF line on
    tabs, whereas mity normalise does

    :param new_vcf: new_vcf is a list of strings, created by merge
    :param out_file: the resulting filename. this should end in vcf.gz
    :return: None. This function writes a vcf.gz and vcf.gz.tbi file.
    """
    f = tmp_mity_file_name()
    logging.debug("Writing uncompressed vcf to " + f)
    with open(f, mode='wt', encoding='utf-8') as myfile:
        for vcf_line in new_vcf:
            myfile.write(vcf_line + '\n')
    gsort_vcf(f, out_file, genome_file=genome_file)

def create_genome_file(vcf_file, genome_file):
    """
    gsort (https://github.com/brentp/gsort) requires a '.genome'
    file to tell it how to sort the vcf records. This function creates a
    '.genome' file in the same order as the contig lines in the vcf header.

    :param vcf_file: a vcf file with the correct contig names
    :param genome_file: the resulting .genome file
    :return: None. this creates a '.genome' file
    """

    vcf = pyvcf.Reader(filename=vcf_file)
    with open(genome_file, mode='wt', encoding='utf-8') as genome_file:
        for contig in vcf.contigs.keys():
            genome_file.write('\t'.join(vcf.contigs[contig]))

def check_dependency(dep, exit=True):
    """
    Check if a dependency exists.

    Special cases:
    * gsort: this utility from brentp has the same name as GNU sort on some linux
      systems. This function will check that brentp's gsort is found.

    >>> check_dependency("ls")
    >>> check_dependency("freebayes")
    >>> check_dependency("gsort")
    :param dep: name of the dependency
    :param exit: If True, then if the dependency isn't found, the session will exit.
    :return: True/False if dependency was found.
    """
    found = False
    try:
        res = subprocess.run(['which', dep], capture_output=True, check=True, encoding="UTF8")
        found = True
        logging.info("Found dependency: " + str(res.stdout).strip())
        if dep == "gsort":
            # There is a potential name clash with gnu sort on linux
            res = subprocess.check_output('gsort --help | grep -c GENOME', shell=True)
            if res != b'2\n':
                logging.error("Adjust your PATH to ensure that brentp's gsort is found before GNU sort" + dep)
    except subprocess.CalledProcessError:
        logging.error("Missing dependency: " + dep)
        if exit:
            sys.exit(1)
    return found


def check_dependencies(f='verchew.ini'):
    """
    Check that all of the dependencies that are listed in the INI format exist.
    :param f: an INI formatted file. see verchew [https://verchew.readthedocs.io/en/latest/]
    :return:
    """
    config = configparser.ConfigParser()
    config.read(f)
    for section in config.sections():
        dependency = config[section]['cli']
        check_dependency(dependency)


def make_hgvs(pos, ref, alt):
    # make HGVS syntax
    if len(alt) > 1 or len(ref) > 1:
        # this is an indel
        if len(ref) > len(alt):
            # this is a del
            delet = ref[1:]
            if len(delet) == 1:
                hgvs_pos = int(pos) + 1
                # print(hgvs_pos)
            elif len(delet) > 1:
                hvgs_pos_start = int(pos) + 1
                hvgs_pos_end = int(pos) + len(delet)
                hgvs_pos = str(hvgs_pos_start) + "_" + str(hvgs_pos_end)
                # print(hgvs_pos)
            # print(hgvs_pos)
            hgvs = "m." + str(hgvs_pos) + "del"
            # print(hgvs)

        else:
            # print("ins")
            # this is an ins
            ins = alt[1:]
            if len(ins) == 1:
                hgvs_pos = int(pos) + 1
                # print(hgvs_pos)
            elif len(ins) > 1:
                hvgs_pos_start = int(pos) + 1
                hvgs_pos_end = int(pos) + len(ins)
                hgvs_pos = str(hvgs_pos_start) + "_" + str(hvgs_pos_end)
                # print(hgvs_pos)
            # print(hgvs_pos)
            hgvs = "m." + str(hgvs_pos) + "ins"
            # print(hgvs)

    else:
        # this is a SNP
        hgvs = "m." + str(pos) + str(ref) + ">" + str(alt)
    return hgvs

def select_reference_fasta(reference, custom_reference_fa=None):
    """
    Allow the user to select one of the pre-loaded reference genome fasta files, via --reference,
    or supply their own via --custom_reference. This function will return the path to
    the reference genome fasta file.

    :param reference: one of the inbuilt reference genomes: hs37d5, hg19, hg38, mm10.
    :param custom_reference_fa: the path to a custom reference genome, or None. If this
    file exists, then it will override the option provided by 'reference'.
    :return the path to the reference genome as a str.

    >>> select_reference_fasta('hg19', None)
    'reference/hg19.chrM.fa'
    >>> select_reference_fasta('hg19', 'mitylib/reference/hs37d5.MT.fa')
    'reference/hs37d5.MT.fa'
    >>> select_reference_fasta('hg19', 'nonexistent.fa')
    'reference/hg19.chrM.fa'

    """
    if custom_reference_fa is not None and os.path.exists(custom_reference_fa):
        res = custom_reference_fa
    else:
        ref_dir = os.path.join(get_mity_dir(), 'reference')
        res = glob('{}/{}.*.fa'.format(ref_dir, reference))
        logging.debug(",".join(res))
        assert len(res) == 1
        res = res[0]
    return res

def select_reference_genome(reference, custom_reference_genome=None):
    """
    This function returns the path to a '.genome' file [1], which is needed for `gsort` to order the chromosomes
    properly.

    [1]: https://bedtools.readthedocs.io/en/latest/content/general-usage.html#genome-file-format

    :param reference: one of the inbuilt reference genomes: hs37d5, hg19, hg38, mm10
    :param custom_reference_genome: the path to a custom reference .genome file, or None. If this
    file exists, then it will override the option provided by 'reference'.
    :return the path to the reference .genome file as a str.

    >>> select_reference_genome('hg19', None)
    'reference/hg19.genome'
    >>> select_reference_genome('hg19', 'mitylib/reference/hs37d5.MT.fa')
    'reference/hs37d5.genome'
    >>> select_reference_genome('hg19', 'nonexistent.fa')
    'reference/hg19.genome'

    """
    if custom_reference_genome is not None and os.path.exists(custom_reference_genome):
        res = custom_reference_genome
    else:
        ref_dir = os.path.join(get_mity_dir(), 'reference')
        logging.debug("Looking for .genome file in " + ref_dir)
        res = glob('{}/{}.genome'.format(ref_dir, reference))
        logging.debug(",".join(res))
        assert len(res) == 1
        res = res[0]
    return res

def get_mity_dir():
    path = os.path.dirname(sys.modules['mitylib'].__file__)
    return path

def vcf_get_mt_contig(vcf):
    """
    get the mitochondrial contig name and length from a VCF file
    :param vcf: path to a vcf file
    :return: a tuple of contig name as str and length as int

    >>> vcf_get_mt_contig('./151016_FR07959656.dedup.realigned.recalibrated.chrMT.dedup.realigned.recalibrated.chrMT.mity.vcf.gz')
    ('MT', 16569)
    """
    r = pyvcf.Reader(filename=vcf, compressed=True)
    chroms = r.contigs.keys()
    mito_contig = {'MT', 'chrM'}.intersection(chroms)
    assert len(mito_contig) == 1
    mito_contig = ''.join(mito_contig)

    return r.contigs[mito_contig].id, r.contigs[mito_contig].length


def bam_get_mt_contig(bam, as_string=False):
    """
    get the mitochondrial contig name and length from a BAM file
    :param bam: path to a bam or cram file
    :return: a tuple of contig name as str and length as int

    >>> bam_get_mt_contig('NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.chrM.bam', False)
    ('chrM', 16569)
    >>> bam_get_mt_contig('NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.chrM.bam', True)
    'chrM:1-16569'
    """
    r = pysam.AlignmentFile(bam, "rb")
    chroms = [str(record.get("SN")) for record in r.header['SQ']]
    mito_contig = {'MT', 'chrM'}.intersection(chroms)
    assert len(mito_contig) == 1
    mito_contig = ''.join(mito_contig)
    res = None
    for record in r.header['SQ']:
        if mito_contig == record['SN']:
            res = record['SN'], record['LN']
    if res is not None and as_string:
        res = res[0] + ":1-" + str(res[1])
    return res

def bam_has_RG(bam):
    """
    Does the BAM or CRAM File have an @RG header? This is critical for mity to correctly call variants.

    :param bam: str: path to bam or cram file
    :return: True/False
    >>> bam_has_RG('NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.chrM.bam')
    """
    r = pysam.AlignmentFile(bam, "rb")
    return len(r.header['RG']) > 0

def get_annot_file(f):
    #mitylibdir = os.path.dirname(inspect.getfile(mitylib))
    #mitylibdir = os.path.dirname(mitylib.__file__)
    mitylibdir = get_mity_dir()
    p = os.path.join(mitylibdir, "annot", f)
    assert os.path.exists(p)
    return p
