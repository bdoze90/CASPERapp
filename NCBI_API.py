"""Parsing GBFF, GFF, and tabular genome annotations.  GFF and tabular do not have the nucleotide sequence while GBFF
has everything"""

from Bio import SeqIO, Entrez
from bs4 import BeautifulSoup
import requests
from ftplib import FTP
import time
import GlobalSettings
import gzip
from PyQt5 import QtWidgets

Entrez.email = "bmendoz1@vols.utk.edu"


class GBFF_Parse:

    def __init__(self, filename, NCBI=False):
        self.filename = filename

    def convert_to_fasta(self):
        count = SeqIO.convert(self.filename, "genbank", self.filename[:-5]+"fasta", "fasta")


class Assembly:

    def __init__(self):
        self.gca_rectList = list()
        self.database_url_list = list()
        self.orgName_dict = dict()
        #Moved the code from this to its own function so that I can return the database URL

    ####### Gets the genbank .gbff NCBI annotation file
    def getNcbiAnnotation(self, organism, databse, ncbi_ret_max):
        ################## start of CASPERweb NCBI API ###################################
        start = time.time()
        
        userInput = organism
        handle = Entrez.esearch(db="assembly", term=userInput , idtype="acc", retmax=ncbi_ret_max)
        record = Entrez.read(handle)

        organismData = []
        linksList = []

        NumOrganisms = len(record['IdList'])
        if NumOrganisms == 0:
            return -1
  
        #creating a string for all the uids
        uidSearch = record['IdList'][0]
        for uid in record['IdList'][1:]:
            uidSearch = uidSearch + ',' + uid
        # print('uid search: ', uidSearch)

        #grabbing all records for all organisms using their uid
        handle = Entrez.efetch(db="assembly", id=uidSearch, rettype="docsum")
        record1 = Entrez.read(handle)

        # print(record1)

        for x in range(NumOrganisms):
            # print("------------------------- NEW RECORD ----------------------")
            # print(record1['DocumentSummarySet']['DocumentSummary'][x])
            if record1['DocumentSummarySet']['DocumentSummary'][x]['FtpPath_GenBank']:
                myDict = {}
                myDict['Organism'] = record1['DocumentSummarySet']['DocumentSummary'][x]['Organism']
                myDict['AssemblyAccession'] = record1['DocumentSummarySet']['DocumentSummary'][x]['AssemblyAccession']
                myDict['genbank_link'] = record1['DocumentSummarySet']['DocumentSummary'][x]['FtpPath_GenBank']
                organismData.append(myDict)
            else:
                continue


        end = time.time()
        

        print("new ncbi api time: ", end - start)
        return (organismData)

    #this function gets a list of GCA's and a list of database URLs for downloading .gff and feature table files
    def getDataBaseURL(self, organism, database, ncbi_ret_max):
        start = time.time()
        # search entrez
        error = False
        print("searching ncbi for: ", organism)

        # this searches entrez searcher. We always search the assembly data base
        # ret-max is taken from the user, my code defaults it to 20 elsewhere
        # and the term is the search organism
        handle = Entrez.esearch(db="assembly", retmax=ncbi_ret_max, term=organism)
        # this parses the data from the search
        record = Entrez.read(handle)
        self.database = database

        # make sure to clear all the things
        self.database_url_list = list()
        self.orgName_dict = dict()
        self.gca_rectList = list()
        self.orgIDs = list()
        myidlist = list()

        # get the internal ID's
        # record is a big dictionary. 'IdList' is a key in that
        for id in record["IdList"]:
            myidlist.append(id)

        # if the len of myIdList is still 0, then return out
        if len(myidlist) == 0:
            return (self.database_url_list, self.orgName_dict)

        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(30)
        # go through and get the GCA/GCF ID's
        for ret in myidlist:
            # this calls Entrez function
            handle = Entrez.esummary(db="assembly", id=ret)
            record = Entrez.read(handle, validate=False)

            # get the orgID which is the ID for the genbank or Refseq link. This could be different than the accession link
            if database == 'RefSeq':
                orgID = record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['RefSeq']
            elif database == 'GenBank':
                orgID = record['DocumentSummarySet']['DocumentSummary'][0]['Synonym']['Genbank']

            # get the accession link and store it
            gca_rec = record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            self.gca_rectList.append(gca_rec)
            self.orgIDs.append(orgID)
            handle.close()
            # sleep so NCBI doesn't kick us out
            time.sleep(0.5)

        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(50)
        # for each GCA_ID, go through and get the refseq/genbank link, and the organism name
        # this is the beautiful soup part
        for i in range(len(self.gca_rectList)):
            url = "https://www.ncbi.nlm.nih.gov/assembly/" + self.gca_rectList[i] + "/"
            source = requests.get(url)
            plain_text = source.text
            soup = BeautifulSoup(plain_text, "html.parser")

            refseq_link = str(soup.find('a', text="FTP directory for RefSeq assembly"))
            genbank_link = str(soup.find('a', text="FTP directory for GenBank assembly"))
            orgName = soup.find('dd')
            refseq_link = refseq_link[refseq_link.find("=") + 2: refseq_link.find(">") - 1]
            genbank_link = genbank_link[genbank_link.find("=") + 2: genbank_link.find(">") - 1]


            if self.database == "GenBank":
                # check and see if GCF is in the the GCA_ID, if so, swap it with GCA
                # not sure if doing this is correct or not, but this way each description actually goes to a download link
                if "GCF" in self.gca_rectList[i]:
                    self.gca_rectList[i] = self.gca_rectList[i].replace("GCF", "GCA")
                database_url = genbank_link
            else:
                database_url = refseq_link

            # check the links and catch errors
            if (database == "RefSeq" and len(refseq_link) < 5):
                error = True
                print("Error: No RefSeq file to download!")
                print(url)
            elif (database == "GenBank" and len(genbank_link) < 5):
                error = True
                print("Error: No GenBank file to download!")
                print(url)
            elif (len(genbank_link) < 5 and len(refseq_link) < 5):
                error = True
                print("Error: No RefSeq or GenBank files to download")
                print(url)
            else:
                error = False
                self.database_url_list.append(database_url)
            # only set the data if there actually is a link to download with
            if orgName and error == False:
                self.orgName_dict[orgName.string + "::" + self.orgIDs[i]] = self.orgIDs[i]

        GlobalSettings.mainWindow.ncbi_search_dialog.searchProgressBar.setValue(80)
        end = time.time()
        
        # print("PART1 TIME: ", end - start)
        # print("database url list: ", self.database_url_list)
        # print("orgname dict: ", self.orgName_dict)
        
        return (self.database_url_list, self.orgName_dict)


    # this function downloads the .gz file from the given link, and stores it in the CSPR_DB folder from Global Settings
    # NOTE: this function should only be used to download .gz files, nothing else
    # NOTE2: this function also only downloads the fna files
    def download_compressed_file(self, database_link):
        print("Downloading: ", database_link)
        ftpLink = database_link[27:]
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        ftp.cwd(ftpLink)

        fileNameIndex = database_link.rfind('/') + 1
        filename = database_link[fileNameIndex:]
        filename = filename + "_genomic.fna.gz"

        local_filename = GlobalSettings.CSPR_DB + "/" +  filename #where the file is stored on the computer
        lf = open(local_filename, "wb")

        ftp.retrbinary("RETR " + filename, lf.write, 8 * 1024)
        lf.close()
        ftp.close()
        print("\tDone downloading")
        return local_filename

    # this function downloads the .gz file from the given link, and stores it in CSPR_DB folder from Global Settings
    # Note: this function should only be used to download .gz files and nothing else
    # Note2: this function will only be able to download feature_tables, gff, gbff files only.
    def download_compressed_annotation_file(self, database_link, fileType):
        # get to the folder on the ftp link
        print("Downloading: ", database_link)
        ftpLink = database_link[27:]
        ftp = FTP('ftp.ncbi.nlm.nih.gov')
        ftp.login()
        ftp.cwd(ftpLink)

        # get the right file
        fileNameIndex = database_link.rfind('/') + 1
        filename = database_link[fileNameIndex:]
        filename = filename + fileType

        # open the right file in the correct location
        local_filename = GlobalSettings.CSPR_DB + "/" + filename
        lf = open(local_filename, "wb")

        # try to write the file, if it is not possible, throw an error
        # write it, and close everything up. Then return the file's location
        try:
            ftp.retrbinary("RETR " + filename, lf.write, 8 * 1024)
        except:
            lf.close()
            ftp.close()
            return
        lf.close()
        ftp.close()
        print("\tDone Downloading")
        return local_filename

    # this function decompresses .gz files into .fna files
    def decompress_file(self, file_path):
        # get the file name
        storeFileName = file_path[:-3]
        print("Decompressing: ", file_path)

        # open the file streams
        readStream = gzip.open(file_path, 'rb')
        writeStream = open(storeFileName, 'w')

        # read the compressed data and close that file
        zipFileContent = readStream.read()
        readStream.close()

        # write that data, and close that file
        tempString = str(zipFileContent)[2:-1]
        tempString = tempString.replace("\\t", "\t")
        tempList = tempString.split("\\n")
        for i in range(len(tempList)):
            writeStream.write(str(tempList[i]) + "\n")
        writeStream.close()

        print("\tDone decompressing")
        return storeFileName

    # this function decompressed an annotation file based on the type of annotation file given
    # stores the file in the CSPR_DB folder, and keeps the same name, just minus the .gz
    # may not need to check for the type of file, not sure. Still need to test
    def decompress_annotation_file(self, file_path, file_type):
        # get the compressed data
        storeFileName = file_path[:-3]
        print("Decompressing: ", file_path)
        readStream = gzip.open(file_path, 'rb')
        writeStream = open(storeFileName, 'w')
        zipFileContent = readStream.read()
        readStream.close()

        # check the file type, and decompress accordingly
        # i may be able to take the if statements out, and decompresss all of them in the same way.
        # just need to talk to brian about that
        if file_type == "_genomic.gbff.gz":
            tempString = str(zipFileContent)[2:-1]
            tempString = tempString.replace("\\t", "\t")
            tempList = tempString.split("\\n")
            for i in range(len(tempList)):
                writeStream.write(str(tempList[i] + "\n"))
        elif file_type == "_genomic.gff.gz":
            tempString = str(zipFileContent)[2:-1]
            tempString = tempString.replace("\\t", "\t")
            tempList = tempString.split("\\n")
            for i in range(len(tempList)):
                writeStream.write(str(tempList[i] + "\n"))
        elif file_type == "_feature_table.txt.gz":
            tempString = str(zipFileContent)[2:-1]
            tempString = tempString.replace("\\t", "\t")
            tempList = tempString.split("\\n")
            for i in range(len(tempList)):
                writeStream.write(str(tempList[i] + "\n"))

        writeStream.close()
        print("\tDone Decompressing")
        return(storeFileName)
# testing below
#myNCBI = Assembly()
#myNCBI.get_annotation_file("Bacillus Subtilis[Organism]", "GenBank")
