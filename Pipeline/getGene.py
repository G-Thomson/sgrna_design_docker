#! /usr/bin/env python

import warnings
import sys
import gffutils
from pyfaidx import Fasta

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def getGene(transcript, genome = "C.sinensis_v2.0_HZAU_csi.chromosome.fa", database = "citrus.db", f=sys.stdout, mode='CDS'):
    
    db = gffutils.FeatureDB(database)
    
    transcript = db[transcript]
    
    if(mode == 'CDS'): # Will return a SeqRecord
        seq=""
        
        if(transcript.strand == '+'):
            for cds in db.children(transcript, featuretype='CDS', order_by='start'):
                seq = seq + cds.sequence(genome, use_strand=True)
            record = SeqRecord(Seq(seq), id=transcript.id, description="")
            return record
        elif(transcript.strand == '-'):
            for cds in db.children(transcript, featuretype='CDS', order_by='start', reverse=True):
                seq = seq + cds.sequence(genome, use_strand=True)
            record = SeqRecord(Seq(seq), id=transcript.id, description="")
            return record
        else:
            warnings.warn("Transcript does not have a strand", UserWarning)
            return
    
    elif(mode == "exons"): # Will return a SeqRecord

        ex_list = []

        if(transcript.strand == '+'):
            for cds in db.children(transcript, featuretype='CDS', order_by='start'):
                seq = cds.sequence(genome, use_strand=True)
                record = SeqRecord(Seq(seq), id=transcript.id, description="")
                ex_list.append(record)
            return ex_list
        elif(transcript.strand == '-'):
            for cds in db.children(transcript, featuretype='CDS', order_by='start', reverse=True):
                seq = cds.sequence(genome, use_strand=True)
                record = SeqRecord(Seq(seq), id=transcript.id, description="")
                ex_list.append(record)
            return ex_list
        else:
            warnings.warn("Transcript does not have a strand", UserWarning)
            return ex_list

# if __name__ == '__main__':

    # a = getGene('orange1.1t03773.1', mode='CDS')
    # print(a.format('fasta'))
    # b = getGene('orange1.1t03773.1', mode='exons')
#     for r in b:
#         print(r.format('fasta'))

    # print(isinstance(a, SeqRecord))
    # print(isinstance(b, list))

    # x = getGene('Cs2g15310.1', mode='CDS')
    # print(x.format('fasta'))
    # y = getGene('Cs2g15310.1', mode='exons')
    # for r in y:
    #     print(r.format('fasta'))
    # from Bio import SeqIO
    # SeqIO.write(y, "example.fasta", "fasta")


# Cs2g15310 rev

# One function
# Create Seq Record
# 
# Return Seq record or list of seq records
# 
# transcript, database, fasta file, output file
# Write function 

# Add annotations
# Add gene_level