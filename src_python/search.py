# Implementing Search in Python

#1-1 Group Name

#1-2 Naive Search Algorithm
import iv2py as iv

def import_fasta_seq_from(filename:str):
    """
    :param filename: name of the file in the 'data' directory (just the filename, not the path)
    :returns: a list of sequences of each record
    """
    records = []
    for record in iv.fasta.reader(file="../data/"+filename):
        records.append(record.seq)
    return records

def naive_search_single(substring:str, string:str):
    """
    :param substring: the substring to be searched for
    :param string: the reference string in which the substring will be searched
    :returns: list of int, that contains the starting index of all occurrences of the substring in the reference string
    """
    store = []
    while True:
        i = string.find(substring)
        if i==-1:
            break
        store.append(i)
        string=string.replace(substring,"%"*len(substring),1)
    return store

def naive_search(reference_filename:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This should have only one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of list of int containing all the starting indices of each occurrence of the read in the reference. 
    """
    reference_sequence = import_fasta_seq_from(reference_filename)
    assert len(reference_sequence)==1, "Reference file contains more than one record"
    reference_sequence = reference_sequence[0]

    read_sequences = import_fasta_seq_from(reads_filename)

    searches = []
    if queries:
        for i in range(queries):
            read = read_sequences[i]
            searches.append(naive_search_single(read,reference_sequence))
    else:
        for read in read_sequences:
            searches.append(naive_search_single(read,reference_sequence))
    return searches

#1-3 Suffix Array Based Search
def binary_search(substring:str, suffix_array:list[int], reference:str):
    """
    :param substring: substring to be searched in the reference
    :param suffix_array: suffix array of the reference string
    :param reference: reference string in which the substring will be searched
    :returns: the index in the suffix array which contains the index where the substring occurs in the reference
    """
    right = len(suffix_array)
    left = 0
    while left<right:
        middle = (right+left)//2
        if reference[suffix_array[middle]:] < substring:
            left = middle + 1
        else:
            right = middle
    return left

import copy
def suffix_search_single(substring, suffix_array, reference):
    """
    :param substring: substring to be searched in the reference
    :param suffix_array: suffix array of the reference string
    :param reference: reference string in which the substring will be searched
    :returns: all the indices where the substring occurs in the reference
    """
    first = binary_search(substring,suffix_array,reference)

    left = first
    while left > 0:
        i = suffix_array[left]
        if reference[i: i + len(substring)] != substring:
            left+=1
            break
        else:
            left -=1

    right = first
    while right < len(suffix_array):
        i = suffix_array[right]
        if reference[i: i + len(substring)] != substring:
            break
        else:
            right +=1
    answer = copy.deepcopy(suffix_array[left:right])
    answer.sort()
    return answer

def suffix_search(reference_filename:str="hg38_partial.fasta.gz", reads_filename:str="illumina_reads_100.fasta.gz", queries:int=None):
    """
    :param reference_filename: name of the file in the 'data' directory that stores the reference sequence. This should have only one record in the file.
    :param reads_filename: name of the file in the 'data' directory that stores the reads sequences. This may have very many records.
    :param queries: if entered, only the first n reads will be searched.
    :returns: list of int containing the starting indices of each read in the reference. -1 if not found.
    """
    reference_sequence = import_fasta_seq_from(reference_filename)
    assert len(reference_sequence)==1, "Reference file contains more than one record"
    reference_sequence = reference_sequence[0]
    suffix_array = iv.create_suffixarray(reference_sequence)

    read_sequences = import_fasta_seq_from(reads_filename)

    searches = []
    if queries:
        for i in range(queries):
            read = read_sequences[i]    
            searches.append(suffix_search_single(read,suffix_array,reference_sequence))
    else:
        for read in read_sequences:
            searches.append(suffix_search_single(read,suffix_array,reference_sequence))
    return searches


#1-4 Benchmark (runtime and memory) your solutions for 1’000, 10’000, 100’000 1’000’000 queries of length 100.

#1-5 Benchmark (runtime) queries of the length 40, 60, 80, and 100 with a suitable number of queries.
