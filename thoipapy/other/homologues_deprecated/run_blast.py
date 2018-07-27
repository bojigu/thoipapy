from Bio.Blast import NCBIWWW

def extract_NCBI_homolgous_as_xml(acc,s,logging,pathdict):
    filename = "OLD script"
    protein=filename[-12:-6]
    #for uniprot_acc in dfset['Protein']:
    #if(uniprot_acc == protein):
    #sys.stdout.write(uniprot_acc)
    #index1=dfset['Protein'].tolist().index(protein)
    #TMD_start=dfset['TMD_start'][index1]
    #TMD_end=dfset['TMD_end'][index1]
    #if protein == "P02724":
    myEntrezQuery = "Homo sapiens[Organism]"
    fasta_string = open(filename).read()
    result_handle = NCBIWWW.qblast("blastp", "nr", fasta_string, expect=10,  hitlist_size=10000)
    myblast="/scratch2/zeng/pythonscript/%s.xml" %protein
    save_file = open(myblast, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    E_VALUE_THRESH = 0.02
    #xml_result_handle = open(myblast)
    #blast_record = NCBIXML.read(xml_result_handle)
    #blast_record = NCBIXML.parse(xml_result_handle)
    #i=0
    #for description in blast_record.descriptions:
        #sys.stdout.write(description.score)
    '''
    #if protein=="P02724":
        #for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < E_VALUE_THRESH:
                    i=i+1
                    query_start=hsp.query_start
                    query_end=hsp.query_end
                    query_seq=hsp.query
                    sbjt_seq=hsp.sbjct
                    tmstr=''
                    sbjstr=''
                    if query_start<=TMD_start and query_end >= TMD_end:
                        sys.stdout.write('****Alignment****',i)
                        #sys.stdout.write('>', alignment.title)
                        #sys.stdout.write('length:', alignment.length)
                        #sys.stdout.write('identity:', hsp.identities)
                        #sys.stdout.write('query_start:', hsp.query_start)
                        #sys.stdout.write('query_end:', hsp.query_end)
                        #sys.stdout.write('subject_start:', hsp.sbjct_start)
                        #sys.stdout.write('subject_end:', hsp.sbjct_end)
                        #sys.stdout.write('e value:', hsp.expect)
                        #sys.stdout.write(hsp.query)
                        #sys.stdout.write(hsp.match)
                        #sys.stdout.write(hsp.sbjct)
                        tm_str_start=TMD_start-query_start
                        tm_str_end=TMD_end-query_start+1
                        k=0
                        j=0
                        tmstr=""
                        sbjtstr=""
                        for char in query_seq:
                            if char != '-':
                                if j >= tm_str_start and j <tm_str_end:
                                    tmstr+=query_seq[j]
                                    sbjtstr+=sbjt_seq[k]
                                j=j+1
                            k=k+1
                        sys.stdout.write(tmstr)
                        sys.stdout.write(sbjtstr)'''