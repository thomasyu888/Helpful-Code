import pysam
import os
import time
import numpy as np
import math
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import itertools
import collections
import pandas
import logging
import pybedtools
import pickle
import tempfile
import time
import re
import yaml
import decimal
from decimal import Decimal

def get_mappedFile_prefix(mapped_file,mode='r',output_dir=None):
    """
    get prefix for a sam/bam file
    output_dir : is used to append that dir as a dir name behind the file prefix
               : useful to direct the bam/sam file to a diff direc
    """
    SUB='get_mappedFile_prefix'
    
    if output_dir:
        dirname = output_dir
    else:
        dirname = os.path.dirname(mapped_file)
    if dirname == "": ##fix when script is run from the same dir where the file is
        dirname = '.'
    if mapped_file.endswith('.bam'):
        file_prefix = dirname + '/' +  os.path.basename(mapped_file).replace('.bam','')
        return (file_prefix)
    elif mapped_file.endswith('.sam'):
        file_prefix = dirname + '/' + os.path.basename(mapped_file).replace('.sam','')
        return (file_prefix)
    else:
        print '[%s]: cant understand the file extension..tried looking for .sam  and .bam' % SUB

    
def get_mappedFile_FH(mapped_file,mode='r'):
    """
    get apt file handle for reading sam / bam file
    """
    SUB='get_mappedFile_FH'
    
    if mapped_file.endswith('.bam'):
        file_handle = pysam.Samfile(mapped_file,('%sb' % mode))
        return (file_handle)
    elif mapped_file.endswith('.sam'):
        file_handle = pysam.Samfile(mapped_file,('%s' % mode))
        return (file_handle)
    else:
        print '[%s]: cant understand the file extension..tried looking for .sam  and .bam' % SUB

def coordinate_sort_bam(bamFile,parallel=False, threads = 2,out_q=None):
    """
    Taking a bamFile as input
    produce a coordinate sorted bam file
    """
    SUB = 'coordinate_sort_bam'
    bamFile_coordinate_sorted_prefix = bamFile.replace('.bam','_cord_sorted')
    bamFile_coordinate_sorted        = bamFile_coordinate_sorted_prefix + '.bam'
    if os.path.exists(bamFile_coordinate_sorted):
        print '[%s]: Bam file %s already coordinate sorted ' % (SUB,bamFile)
        return bamFile_coordinate_sorted
    else:
        print '[%s]: Name sorting bam file %s' % (SUB,bamFile)
        if parallel:
            args = ['sambamba','sort','-t',str(threads),'-o',bamFile_coordinate_sorted,bamFile]
        else:
            args = ['samtools','sort',bamFile,bamFile_coordinate_sorted_prefix]

        return_code = subprocess.check_call(args)  
        if return_code == 0:
            print '[%s]: Created coordinate sorted bam for %s' % (SUB,os.path.basename(bamFile))
            if out_q:
                out_q.put(bamFile_coordinate_sorted)
            return bamFile_coordinate_sorted
        else:
            print '[%s]: Error creating coordinate sorted bam for %s' % (SUB,os.path.basename(bamFile))
            if out_q:
                out_q.put(False)

  
def bam_to_bedGraph(bamFile, scale=False, sorted=False, **kwargs):
    SUB='bam_to_bedGraph'
    if sorted is not True:
        #sort the bam file by coordinates
        bamFile_sorted = coordinate_sort_bam(bamFile)

    base_path = os.path.dirname(bamFile) or '.'
    bamFile_name = os.path.basename(bamFile)
    bam_prefix = get_mappedFile_prefix(bamFile)
    out_bedgraphFile = base_path + '/' + bamFile_name.replace('.bam','') + '.bedgraph'
    
    if os.path.exists(out_bedgraphFile):
        print '[%s]: Expected output %s already present no processing required' % (SUB,out_bedgraphFile)
        return(out_bedgraphFile)

    #scale the read counts if asked
    if scale is True:
        readcount = get_mapped_read_count(bamFile_sorted)
        factor = 1 / (readcount / 1e6)
    else:
        factor = 1
     
    #bedgraph option line
    bedgraph_line = 'track type=bedGraph name=%s color=%s altColor=%s' % (bam_prefix,'43,131,186','171,121,164') 
    #find the genome coverage
    bam_bedtool = pybedtools.BedTool(bamFile_sorted)
    bam_bedtool.genome_coverage(bga=True, scale=factor,trackopts=bedgraph_line, output=out_bedgraphFile)
    
    print '[%s]: Converted bam %s to bedgraph file'  % (SUB,bamFile)
    return out_bedgraphFile


def bam_to_bigWig(bamFile,scale=False, sorted=False, **kwargs):
    """
    For a given bam file
    create a bigWig File
    """
    SUB = "bam_to_bigWig"
    #user supplied args if any
    force = kwargs.get('force',False)
    output_dir = kwargs.get('output_dir',False)

    #output bigWig filename
    out_bigWig_file = get_mappedFile_prefix(bamFile,output_dir=output_dir)  + '.bw'    
    if os.path.exists(out_bigWig_file):
        print '[%s]: Expected output %s already present no processing required' % (SUB,out_bigWig_file)
        return(out_bigWig_file)

    #constructing genome size file on the fly from bam file
    bam_fh = get_mappedFile_FH(bamFile)
    bam_prefix = get_mappedFile_prefix(bamFile)
    temp_genome_file = tempfile.NamedTemporaryFile(mode='w+t',suffix='.tmp', delete=False) #when delete is True file is deleted as soon as it closed
    [ temp_genome_file.write('%s\t%s\n' % (chr,len)) for chr,len in itertools.izip(bam_fh.references,bam_fh.lengths) ]
    temp_genome_file.close()


    #1. convert bam to bedgraph first 
    bedGraph_file = bam_to_bedGraph(bamFile, scale=scale, sorted=False, **kwargs )

    #2. convert to bigwig using the UCSC script   #should be in the path
    cmds = ['bedGraphToBigWig',bedGraph_file, str(temp_genome_file.name),out_bigWig_file]
    p = subprocess.Popen(cmds,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    stdout,stderr = p.communicate()
    
    if not stderr:
        print 'Created a bigWig file %s for %s' % (os.path.basename(out_bigWig_file),
                                                   os.path.basename(bamFile))
    else:
        print stderr
    return(out_bigWig_file)

bams = os.listdir('../combined')
for i in bams:
    if "bam" in i:
        bamFile = "../combined/"+i
        bam_to_bigWig(bamFile,output_dir = ".")
