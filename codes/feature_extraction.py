"""
    Extract the features used in the manuscript
"""


from tqdm import tqdm
import numpy as np
from sklearn.decomposition import PCA
import re
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from math import log2


def extractFeatures(genome_sequence, model_params):
    """
    Extracts the features from the genome sequence
    
    Arguments:
        genome_sequence {string} -- input genome sequence
        model_params {dictionary} -- SSG-LUGIA model configuration
    
    Returns:
        2d numpy.array -- computed features of the windowed sequences
    """

    w = model_params['w']                   # length of sliding window
    dw = model_params['dw']                 # step size of sliding window
    
    X = []                                  # computed features


    dinucleotide_freq = []                  # dinucleotide frequency
    all_aa_freq = []                        # amino acid frequency
    kmer4_freq = []                         # 4-mer frequency    


    kmer2_map = generateKmerList(2)       # 2-mer mapping
    kmer3_map = generateKmerList(3)       # 3-mer mapping
    kmer4_map = generateKmerList(4)       # 4-mer mapping


    amino_acid_mapper = generateAminoAcidMapping(kmer3_map)  # amino acid mapping
    

    cur_kmer2_freq = [0] * len(kmer2_map)       # initialize freuency of 2-mers
    cur_kmer3_freq = [0] * len(kmer3_map)       # initialize freuency of 3-mers
    cur_kmer4_freq = [0] * len(kmer4_map)       # initialize freuency of 4-mers


    lastseq = None              # initializing the last sequence for efficient computation

    print('Extracting Features')
    for st in tqdm(range(0,len(genome_sequence),dw)):

        if((st+w)>=len(genome_sequence)):           
            break

        X.append([])

        window_seq = genome_sequence[st:st+w]           # current window

        cur_kmer2_freq = computeKmerFreq(window_seq, lastseq, dw, 2, cur_kmer2_freq, kmer2_map)   # compute frequency of 2-mers
        cur_kmer3_freq = computeKmerFreq(window_seq, lastseq, dw, 3, cur_kmer3_freq, kmer3_map)   # compute frequency of 3-mers
        cur_kmer4_freq = computeKmerFreq(window_seq, lastseq, dw, 4, cur_kmer4_freq, kmer4_map)   # compute frequency of 4-mers
        
        kmer4_freq.append(np.array(cur_kmer4_freq))

                                                            # computing nucleotide frequency
        nt_cnt = {'A':window_seq.count('A'),
                  'T':window_seq.count('T'),
                  'C':window_seq.count('C'),
                  'G':window_seq.count('G')}
        
        gc_content = computeGCContent(nt_cnt['G'], nt_cnt['C'], w)        # computing gc content
        gc_skew = computeGCSkew(nt_cnt['G'], nt_cnt['C'])                 # computing gc skew

        X[-1].append(gc_content)
        X[-1].append(gc_skew)

        karlin_di = computeKarlinDinucleotide(nt_cnt, cur_kmer2_freq, kmer2_map,mode=model_params['karlin_mode'])         # computing karlin bias

        aa_freq = computeAaFreq(cur_kmer3_freq, kmer3_map, amino_acid_mapper)                 # computing amino acid frequency

        if(model_params['entropy_features']):
                                                                            # computing entropy features
            karlin_di_entropy = entropyFromList(np.array(karlin_di))
            aa_freq_entropy = entropyFromList(np.array(aa_freq))            
            kmer4_freq_entropy = entropyFromList(np.array(cur_kmer4_freq))

            X[-1].append(karlin_di_entropy)
            X[-1].append(aa_freq_entropy)
            X[-1].append(kmer4_freq_entropy)


        dinucleotide_freq.append(karlin_di)
        all_aa_freq.append(aa_freq)


        lastseq = window_seq[:]             # update the last sequence for the next iteration


    if(model_params['pca_dn']>0):                                                   # compute pca components of dinucleotide bias
        pca0 = PCA(n_components=model_params['pca_dn'], svd_solver='full')
        XK_dn = pca0.fit_transform(dinucleotide_freq)

    pca2 = PCA(n_components=model_params['pca_amino_acid'], svd_solver='full')      # compute pca components of amino acid frequency
    XK_aa = pca2.fit_transform(all_aa_freq)

    pca1 = PCA(n_components=model_params['pca_kmer4'], svd_solver='full')           # compute pca components of 4-mer frequency
    XK_4 = pca1.fit_transform(kmer4_freq)
    

    for i in range(len(X)):
        if(model_params['pca_dn']>0):
            X[i].extend(XK_dn[i])
        X[i].extend(XK_aa[i])
        X[i].extend(XK_4[i])        
        
    
    X = np.array(X)
    
    
    return X


def entropyFromList(inp_list):
    """
    Computes entropy from a list
    
    Arguments:
        inp_list {list} -- input list containing frequencies
    
    Returns:
        float -- computed entropy
    """

    inp_arr = np.array(inp_list)        # converting to numpy array for efficient computation
    summ = np.sum(inp_list)

    if(abs(summ)<1e-6):                 # if number too small log will cause error
        return 0

    p_arr = (inp_arr / summ) + 1        # normalizing the array, added 1 to avoid zeroes

    log2_p_arr = np.log2(p_arr)

    entropy =  (np.sum(p_arr * log2_p_arr))        # computing array using the standard formula

    return entropy



def generateKmerList(k):
    """
    Generates a mapping of all the k-mers for a specific k
    
    Arguments:
        k {int} -- 'k' of k-mer
    
    Returns:
        dictionary -- mapping of the kmers against an arbitrary index
    """

    kmers = ['']                            # initializing list of k-mers
    alphabet = ['A', 'T', 'C', 'G']         # list of possible alphabets
    currlen = 0                             # index variable

    while(currlen<k):

        t_kmer = []                         # temporary
        
        for kmer in kmers:
            
            t_kmer.append(kmer+alphabet[0])
            t_kmer.append(kmer+alphabet[1])
            t_kmer.append(kmer+alphabet[2])
            t_kmer.append(kmer+alphabet[3])

        kmers = np.array(t_kmer)

        currlen += 1

    kmer_list = kmers[:]                # list of all the k-mers
    kmer_map = {}                       # initializing the mapping

    for i in range(len(kmer_list)):     

        kmer_map[kmer_list[i]] = i      # setting the index for the k-mer

    return kmer_map

def generateAminoAcidMapping(kmer3_map):
    """
    Generates a mapping of the amino acids
    
    Arguments:
        kmer3_map {dictionary} -- mapping of the 3-mers

    Returns:
        dictionary -- mapping of the amino acids against an arbitrary index
    """

    amino_acid_mapper = {}
    aa_index = 0                                    # initialize index

    for codon in kmer3_map:
        
        coding_dna = Seq(codon, IUPAC.unambiguous_dna)      
        aa = coding_dna.translate(to_stop=False)            # translate codon into amino acid
        aa= str(aa)
        
        amino_acid_mapper[coding_dna] = aa

        if aa in amino_acid_mapper:

            amino_acid_mapper[aa].append(codon)

        else:

            amino_acid_mapper[aa] = [aa_index, codon]
            aa_index += 1

    return amino_acid_mapper


def computeGCContent(cnt_g, cnt_c, lenn):
    """
    Computes GC Content
    
    Arguments:
        cnt_g {int} -- frequency Count of G
        cnt_c {int} -- frequency Count of C
        lenn {int} -- length of genome
    
    Returns:
        float -- computed GC Content
    """
    
    return (cnt_g + cnt_c) / (lenn)         # computing GC Content using the standard formula



def computeGCSkew(cnt_g, cnt_c):
    """
    Computes GC Skew
    For more reference please refer to https://www.ncbi.nlm.nih.gov/pmc/articles/PMC275568/
    
    Arguments:
        cnt_g {int} -- frequency Count of G
        cnt_c {int} -- frequency Count of C
    
    Returns:
        float -- computed GC Skew
    """
    
    return (cnt_g - cnt_c) / (cnt_g + cnt_c)    # computing GC Skew using the standard formula



def computeKarlinDinucleotide(frequency, kmer2_freq, kmer2_map, mode='normalized'):
    """
    Computes Karlin Dinucleotide Bias
    For more reference please refer to https://www.sciencedirect.com/science/article/abs/pii/S1369527498800957

    In addition we have proposed a normalized version of Karlin Dinucleotide Bias
    
    Arguments:
        frequency {dictionary} -- frequency of the nucleotides
        kmer2_freq {list} -- frequency of the 2-mers
        kmer2_map {dictionary} -- mapping of the 2-mers against indices
    
    Keyword Arguments:
        mode {str} -- mode of karlin dinucleotide bias computation (default: {'normalized'})
    
    Returns:
        numpy.array -- computed karlin dinucleotide bias
    """
    
    if(mode=='original'):                       # original formulation of karlin dinucleotide bias

        Karlin_dinucleotide = np.array(kmer2_freq)
        Karlin_dinucleotide = Karlin_dinucleotide.astype(float)         # IMPERATIVE type casting

        for a in ['A','T','C','G']:

            for b in ['A','T','C','G']:
                
                Karlin_dinucleotide[kmer2_map[a+b]] = Karlin_dinucleotide[kmer2_map[a+b]] / (frequency[a] * frequency[b] + 1e-6)    # original formula
                                
        return Karlin_dinucleotide


    elif(mode=='raw'):                          # no computation, merely returns the 2-mer frequency

        Karlin_dinucleotide = np.array(kmer2_freq)
        Karlin_dinucleotide = Karlin_dinucleotide.astype(float)         # IMPERATIVE type casting
        return Karlin_dinucleotide


    elif(mode=='normalized'):                   # proposed normalized formulation of karlin dinucleotide bias

        Karlin_dinucleotide = np.array(kmer2_freq)
        Karlin_dinucleotide = Karlin_dinucleotide.astype(float)         # IMPERATIVE type casting
        
        for a in ['A','T','C','G']:

            for b in ['A','T','C','G']:
                                
                Karlin_dinucleotide[kmer2_map[a+b]] = Karlin_dinucleotide[kmer2_map[a+b]] /  ( ( (frequency[a])**0.5 ) * ( (frequency[b])**0.5 ) + 1e-6) # proposed formula
                                
        return Karlin_dinucleotide



def computeKmerFreq(seq, lastseq, step, kmerlen, kmer_freq, kmer_map):
    """
    Computes the k-mer frequency in a sliding window manner
    Using a deque like idea improves time complexity from O(nk) to O(n)
    
    Arguments:
        seq {string} -- current sequence under consideration
        lastseq {string} -- the immediatly previous sequence, used to update k-mer frequence efficiently
        step {int} -- length of step of sliding window
        kmerlen {int} -- length of 'k' in k-mer
        kmer_freq {list} -- frequency of the k-mers
        kmer_map {dictionary} -- mapping of the k-mers against indices
    
    Returns:
        numpy.array -- computed frequency of the k-mers
    """
        
    if(type(lastseq)==type(None)):              # for first computation, there is no last sequence
                                                # so the whole computation is performed
        for i in range(len(seq)-kmerlen+1):         
            try:
                kmer_freq[kmer_map[seq[i:i+kmerlen]]] += 1                      # update the kmer frequency
            except:             # in case there are some error in sequence
                pass

        
    else:                                       # already computed the frequency for certain part of the current window
                                                # so merely remove the part from the last one and add the part of the current one
        for i in range(step):
            try:
                kmer_freq[kmer_map[lastseq[i:i+kmerlen]]] -= 1          # removing the k-mer from the previous window
            except:
                pass

        
        for i in range(len(seq)-step-kmerlen+1, len(seq)-kmerlen+1):
            try:
                kmer_freq[kmer_map[seq[i:i+kmerlen]]] += 1              # adding the k-mer of the current window
            except:
                pass

    return np.array(kmer_freq)                  # copying into a numpy array



def computeAaFreq(kmer_3_freq, kmer3_map, aa_map):
    """
    Computes the Amino Acid Frequency
    
    Arguments:
        kmer_3_freq {list} -- frequency of the 3-mers
        kmer3_map {dictionary} -- mapping of the 3-mers against indices
        aa_map {dictionary} -- mapping of the amino acids against codons or 3-mers
    
    Returns:
        numpy.array -- computed amino acid frequency
    """

    aa_frq = [0] * 21                       # initialize the amino acid frequency

    for codon in kmer3_map:                 

        aa_frq[aa_map[aa_map[codon]][0]] += kmer_3_freq[kmer3_map[codon]]       # maps the codons into amino acids and updates frequency

    return np.array(aa_frq)                  # copying into a numpy array
