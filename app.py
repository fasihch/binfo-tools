import json
from collections import Counter, OrderedDict

import matplotlib.pyplot as plt
import numpy as np
import requests
import streamlit as st
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils import GC, MeltingTemp, GC_skew, seq3
from streamlit_lottie import st_lottie



def load_lottiefile(filepath: str):
    with open(filepath, "r") as f:
        return json.load(f)


def load_lottieurl(url: str):
    r = requests.get(url)
    if r.status_code != 200:
        return None
    return r.json()


lottie_coding = load_lottiefile("lottie\hello.json")  # replace link to local lottie file
lottie_hello = load_lottieurl("https://assets7.lottiefiles.com/packages/lf20_M9p23l.json")

st_lottie(
    lottie_hello,
    speed=1,
    reverse=False,
    loop=True,
    quality="low",  # medium ; high
    height=None,
    width=None,
    key=None,
)
# Show description and sequence of selected DNA in dot plot section
def show_des_and_seq_for_dot_plot(id1_des, id2_des, id1_seq, id2_seq, unique_key):
    file_details = st.radio("Details", ("Description", "Sequence"), key=unique_key)
    if file_details == "Description":
        st.write("File 1 Description:")
        st.write(id1_des)
        st.write("File 2 Description")
        st.write(id2_des)

    elif file_details == "Sequence":
        st.write("File 1 Sequence")
        st.write(id1_seq)
        st.write("File 2 Sequence")
        st.write(id2_seq)


# Complete all tasks in DNA Analysis section
def complete_tasks(full_seq, des, unique_key):
    file_details = st.radio("Details", ("Description", "Sequence"), key=unique_key)

    # Show description and sequence in DNA Analysis section
    if file_details == "Description":
        st.write(des)
    elif file_details == "Sequence":
        st.write(full_seq)

    # Nucleotide occurances plot and color selector for the bars
    st.subheader("Plot Nucleotide Frequency")
    full_seq_freq = OrderedDict(Counter(full_seq))

    bar1_colour = st.color_picker("Pick Colour for Bar 1", key=unique_key)
    bar2_colour = st.color_picker("Pick Colour for Bar 2", key=unique_key)
    bar3_colour = st.color_picker("Pick Colour for Bar 3", key=unique_key)
    bar4_colour = st.color_picker("Pick Colour for Bar 4", key=unique_key)

    if st.button("Plot Frequency", key=unique_key):
        barlist = plt.bar(full_seq_freq.keys(), full_seq_freq.values())
        barlist[0].set_color(bar1_colour)
        barlist[1].set_color(bar2_colour)
        barlist[2].set_color(bar3_colour)
        barlist[3].set_color(bar4_colour)
        st.pyplot()

    st.subheader("Properties")
    st.set_option('deprecation.showPyplotGlobalUse', False)

    # GC Content, GC Melting temp, GC_skew, Complement and reverse complement
    gc_count = GC(full_seq)
    st.write("GC Content: {}".format(gc_count))

    mt = MeltingTemp.Tm_GC(full_seq, strict=False)
    st.write("Melting Temperature based on GC Content: {}".format(mt))

    gc_skew_bases = st.number_input("Enter number of bases", key=unique_key)
    try:
        gc_skew = GC_skew(full_seq, int(gc_skew_bases))
        st.write("GC Skew for {} bases: {}".format(gc_skew_bases, gc_skew))
    except ValueError:
        st.write("Enter a Valid Number for bases")

    if st.checkbox("Complement", key=unique_key):
        st.write(full_seq.complement())

    elif st.checkbox("Reverse Complement", key=unique_key):
        st.write(full_seq.reverse_complement())

    # Protein Synthesis
    st.subheader("Protein Synthesis")
    p1 = full_seq.translate()
    if st.checkbox("Transcription: DNA to mRNA", key=unique_key):
        st.write(full_seq.transcribe())

    elif st.checkbox("Translation: DNA to 1 letter Amino Acid Sequence", key=unique_key):
        st.write(p1)

    elif st.checkbox("Translation: DNA to 3 letter Amino Acid Sequence", key=unique_key):
        full_aa_name = str(p1).replace("*", "")
        st.write(seq3(full_aa_name))

    elif st.checkbox("Plot Amino Acid Frequency", key=unique_key):
        aa_freq = OrderedDict(Counter(str(p1)))
        bar_colour = st.color_picker("Pick Colour for all Bars", key=unique_key)
        plt.bar(aa_freq.keys(), aa_freq.values(), color=bar_colour)
        st.pyplot()
        st.write("Asterisk (*) - Denotes Stop Codons.")

# delta, matrix, makematrix and dotplotx functions for making the dotplot in the dotplot section
def delta(x, y):
    return 0 if x == y else 1


def matrix(seq1, seq2, i, j, k):
    return sum(delta(x, y) for x, y in zip(seq1[i:i + k], seq2[j:j + k]))


def makeMatrix(seq1, seq2, k):
    n, m = len(seq1), len(seq2)
    return [[matrix(seq1, seq2, i, j, k) for j in range(m - k + 1)] for i in range(n - k + 1)]


def dotplotx(seq1, seq2):
    plt.imshow(np.array(makeMatrix(seq1, seq2, 1)))
    xt = plt.xticks(np.arange(len(list(seq2)), len(list(seq2))))
    yt = plt.yticks(np.arange(len(list(seq1)), len(list(seq1))))
    plt.show()

# Main function where the website is controlled
def main():
    """BInfo Aid"""

    st.title("BInfo AID")

    # To select different pages on the webpage
    menu = ["DNA Analysis", "Dot Plot",]
    choices = st.sidebar.selectbox("Select Activity", menu)

    # Go to DNA Analysis page
    if choices == "DNA Analysis":
        st.subheader("What is BInfo Aid?")
        st.write("It is an easy to use web app for analysing DNA sequences. BInfo aid is a web app that is \
                  designed for extracting the meaningful information from the mass of molecular biology,\
                  biological databases & to carry out sequence or structural analysis.It is a data-mining \
                  software that retrieves data from genomic sequence databases and also visualization tools  \
                  to analyze and retrieve information from proteomic databases.")\


        st.write("Specficially this is a simple webpage to showcase the power of bioinformatics by\
                 displaying basic DNA Analysis that can be done using Python.")

        st.write("In order to use this functionality, download a FASTA file or enter a NCBI ID from the Nucleotide database\
                on the NCBI website.")

        st.subheader("DNA Sequence Analysis")

        # Check for offline upload of files
        if st.checkbox("Offline File Upload"):
            seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])
            if seq_file is not None:
                record = SeqIO.read(seq_file, "fasta")
                complete_tasks(record.seq, record.description, "1")

        st.subheader("Or")

        # Check for Online NCBI ID
        if st.checkbox("Online File"):
            ncbi_id = st.text_input("Enter NCBI ID to Retrieve Data")
            if ncbi_id != "":
                Entrez.email = "fasihahmeduaf@gmail.com"
                handle = Entrez.efetch(db="nucleotide", id=ncbi_id, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                complete_tasks(record.seq, record.description, "2")


    # Go to Dot Plot page
    elif choices == "Dot Plot":
        st.write("In order to use this functionality, download or enter two FASTA files/NCBI IDs from the Nucleotide database\
                on the NCBI website.")

        st.subheader("Dot Plot Generator between 2 Sequences")

        # Check for offline upload of files
        if st.checkbox("Offline File Upload"):
            seq1_file = st.file_uploader("Upload FASTA File 1", type=["fasta", "fa"])
            seq2_file = st.file_uploader("Upload FASTA File 2", type=["fasta", "fa"])

            if (seq1_file != None) and (seq2_file != None):
                record1, record2 = SeqIO.read(seq1_file, "fasta"), SeqIO.read(seq2_file, "fasta")
                full_seq1, full_seq2 = record1.seq, record2.seq

                show_des_and_seq_for_dot_plot(record1, record2, full_seq1, full_seq2, "1")

                user_limit = st.number_input("Enter the First Number of Bases to Compare", 10, 200, 50, key="1")
                st.subheader("Dotplot for first {} Nucleotides".format(user_limit))
                dotplotx(full_seq1[0:user_limit], full_seq2[0:user_limit])
                st.pyplot()

        st.subheader("Or")

        # Check for online NCBI IDs
        if st.checkbox("Online File"):
            ncbi_id1 = st.text_input("Enter NCBI ID 1 to Retrieve Data")
            ncbi_id2 = st.text_input("Enter NCBI ID 2 to Retrieve Data")

            if (ncbi_id1 != "") and (ncbi_id2 != ""):
                Entrez.email = "fasihahmeduaf@gmail.com"
                handle1 = Entrez.efetch(db="nucleotide", id=ncbi_id1, rettype="fasta", retmode="text")
                record1 = SeqIO.read(handle1, "fasta")
                id1_des, id1_seq = record1.description, record1.seq
                handle1.close()

                handle2 = Entrez.efetch(db="nucleotide", id=ncbi_id2, rettype="fasta", retmode="text")
                record2 = SeqIO.read(handle2, "fasta")
                id2_des, id2_seq = record2.description, record2.seq
                handle2.close()

                show_des_and_seq_for_dot_plot(id1_des, id2_des, id1_seq, id2_seq, "2")

                user_limit = st.number_input("Enter the First Number of Bases to Compare", 10, 200, 50, key="2")
                st.subheader("Dotplot for first {} Nucleotides".format(user_limit))
                dotplotx(id1_seq[0:user_limit], id2_seq[0:user_limit])
                st.pyplot()


# Start program when streamlit is run
if __name__ == '__main__':
    main()
