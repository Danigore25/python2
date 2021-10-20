#!/usr/bin/env python

'''NAME

        %(progname)s

VERSION

        %(version)s

AUTHOR

    Daniela Goretti Castillo León <danigore22@gmail.com> <dgoretti@lcg.unam.mx>

DESCRIPTION

    Obtiene el abstract y citas de un artículo de PubMed usando una lista de PMIDs. 

CATEGORY
    Base de datos

USAGE
    %(usage)s
        

ARGUMENTS
    -h, --help            show this help message and exit
    -i #, --input=#       list of pmids # (one identifier by line)
                          if not specified, the standard input is used

EXAMPLES

    Input:   [ -i filename]
        34282943
        32611704
        31406982
        30625167

    Output:
    
        Abstract and citations


SEE ALSO
        [ poner el otro programa relacionado de search author ]

GITHUB LINK
    https://github.com/Danigore25/python2/blob/master/src_homework/%(progname)s
        

'''
# ===========================================================================
# =                            imports
# ===========================================================================
import os
import sys
import argparse
from Bio import Entrez
from Bio import Medline


# ===========================================================================
# =                            functions
# ===========================================================================
def fetch_abstract(idlist):
    '''Fetch the medline record for pmid list

        :param idlist: The pmid list
        :type idlist: string list

        :returns: object with medline records
        :rtype: object list
    '''  
    fetch_handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
    records = list(Medline.parse(fetch_handle))
    fetch_handle.close()

    return records


def get_cited_by(pmid):
    '''Get the pmcid references that the input is cited by

        :param pmid: The pmid identifier
        :type pmid: string

        :returns: pmcid list
        :rtype: list
    ''' 
    get_pmc = Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", from_uid=pmid)
    search_results = Entrez.read(get_pmc)
    if len(search_results[0]['LinkSetDb']) > 0:
        pmc_ids = [link["Id"] for link in search_results[0]["LinkSetDb"][0]["Link"]]
    else:
        pmc_ids = []
    return pmc_ids


def main(args):
    Entrez.email = "dgoretti@lcg.unam.mx"

    # Reading the pmids list from file
    inputfile = open(args.input, "r")
    pmid_list = inputfile.read().splitlines()
    inputfile.close()

    # fetching medline records from pmid list
    pubmed_records = fetch_abstract(pmid_list)

    # Print the fields from medline record list
    for record in pubmed_records:
        id = record.get("PMID", "?")
        print("PMID ",id)
        print("TI   ",record.get("TI", "?"))
        print("AU   ",record.get("AU", "?"))
        print("\n")
        print("AB   ",record.get("AB", "?"))

        # get the cites for each pmid
        cites_list = get_cited_by(id)
        print("CITE ", cites_list)
        print("\n//\n")
        
# ===========================================================================
# =                            main
# ===========================================================================
if __name__ == '__main__':

    USAGE = '''%s -i inputfile -o outputfile [-h]'''
    VERSION = '20211019'
    PROG_NAME = os.path.basename(sys.argv[0])

    # arguments
    doc =  globals()['__doc__'] % {'usage' : USAGE % PROG_NAME, 'version' : VERSION, 'progname' : PROG_NAME}
    parser = argparse.ArgumentParser(usage=USAGE % PROG_NAME, description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input", action="store", dest="input", metavar="#", help="PMIDlist file")
    #parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout, help='output file name')   
    args = parser.parse_args()

    
    if not args.input:
        parser.print_usage()
        sys.exit()
    if len(sys.argv) == 1:
        print(USAGE % PROG_NAME)
        sys.exit(0)    
    try:
        main(args)
    except:
        sys.stderr.write('Error while running\n')
    
