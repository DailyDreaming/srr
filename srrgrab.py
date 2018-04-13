import os
import sys
import subprocess
import time
from subprocess import check_call
from os.path import dirname, realpath

# BioPython
from Bio import Entrez
from Bio.Blast import NCBIWWW # NCBI API wrapper
from Bio.Blast import NCBIXML # parse output BLAST XML
from Bio import SeqIO # ingests fastas


# SRR5383888
# Human w/ HBV (PacBio): https://www.ncbi.nlm.nih.gov/sra/SRX2678953[accn]

# SRR5383891
# Human w/ HBV (PacBio; Control): https://www.ncbi.nlm.nih.gov/sra/SRX2678956[accn]

# SRR6172653
# PacBio Soil Samples: https://www.ncbi.nlm.nih.gov/sra/SRX3283433[accn]

# SRR6172655
# PacBio Soil Samples: https://www.ncbi.nlm.nih.gov/sra/SRX3283431[accn]

def datapath(path):
    return os.path.abspath(os.path.join(os.path.dirname(realpath(__file__)), '..', 'data', path))

def toolpath(path):
    return os.path.abspath(os.path.join(os.path.dirname(realpath(__file__)), '..', 'tools', path))

def homepath(path):
    return os.path.abspath(os.path.join(dirname(os.path.realpath(__file__)), path))

class VirLib(object):
    def __init__(self):
        self.SRA_tool_path = self.installSRAtoolkit()
        self.fastq_dump_path = os.path.join(self.SRA_tool_path, 'fastq-dump')

    def installSRAtoolkit(self):
        """
        Compiled and tested on Ubuntu.

        :return:
        """
        print('Installing SRA Toolkit.')
        tools_dir = os.getcwd()
        SRA_tar_path = os.path.join(tools_dir, 'sratoolkit.current-mac64.tar.gz')
        SRA_tar_url = 'http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz'
        # decompresses into version 2.9.0 for some reason
        SRA_tool_path = os.path.join(tools_dir, 'sratoolkit.current-mac64/bin')
        # /home/lifeisaboutfishtacos/Desktop/virusRepo/LongReadViruses/sratoolkit.2.9.0-centos_linux64/bin

        # if there's no SRA_toolkit tar, fetch it from ftp
        if not os.path.exists(SRA_tar_path):
            ftp_fetch = ["wget", "-P", tools_dir, SRA_tar_url]
            subprocess.check_call(ftp_fetch)
        # extract the compressed SRA_toolkit if not already extracted
        if not os.path.exists(SRA_tool_path):
            decompress_SRA_tools = ['tar', '-xzf', SRA_tar_path]
            subprocess.check_call(decompress_SRA_tools)

        assert os.path.exists(SRA_tool_path)
        return SRA_tool_path

    def getSRR(self, SRR_code):
        """
        Takes an SRR code, downloads it as a fastq, and returns the path to that fastq file.

        API Documentation: https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump

        :param SRR_code:
        :return: filepath of downloaded fasta
        """
        print('Fetching the SRR.')
        # basic example cmd: "./fastq-dump -X 5 -Z SRR5150787"
        sratools_API_cmd = "{fastq_dump_path} -Z {SRR_code} > {SRR_code}.fastq" \
                           "".format(fastq_dump_path=self.fastq_dump_path,
                                     SRR_code=SRR_code)
        # run the command
        print('Now running: {cmd}'.format(cmd=sratools_API_cmd))
        process = subprocess.Popen(sratools_API_cmd, shell=True)
        process.communicate()
        # assert the output file was created and return
        SRR_filepath = os.path.join(os.getcwd(), SRR_code + '.fastq')
        assert os.path.exists(SRR_filepath)

        # throttle requests so that we don't contact the server more than once every 10 seconds
        # More Info: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=DeveloperInfo
        time.sleep(11)

        return os.path.join(SRR_filepath)


v = VirLib()
v.getSRR('SRR5209994')
