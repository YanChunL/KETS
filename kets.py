# This is the code to evaluate simulated genome assembly results based on kmer (include CR, SCR, DCR, APLC, ADF, LKLIS, ICADF)
# -*- coding: UTF-8 -*-

import argparse
import time
from pandas import *
from random import *
import os
import re
from bisect import bisect_left
from collections import defaultdict
def ArgParse():
    group = argparse.ArgumentParser(description='A python script for genome assessment.')
    group.add_argument('-i', '--input', help='assemble result with fasta format.', required=True)
    group.add_argument('-r', '--reference', help='reference sequence with fasta format.', required=True)
    group.add_argument('-k', '--kmer-length', type=int, help='the kmer length used in assessment, default=21.',default=21)
    group.add_argument('-o', '--out-prefix',help='prefix of output files.',required=True)
    group.add_argument('-s', '--sample', help='the number of ref unique kmer sampled, default=all.', default="all")
    return group.parse_args()


def find_gaps(fasta_file):
    with open(fasta_file, 'r') as file:
        sequences = file.read().split('>')
    gap_dic={}
    for seq in sequences[1:]:
        lines = seq.strip().split('\n')
        scaffold_name = lines[0].split()[0]
        sequence = ''.join(lines[1:])
        gaps = [(m.start(), m.end() - 1) for m in re.finditer(r'[nN]+', sequence)]
        gap_dic[scaffold_name]= gaps
    gap_list = [(key, value) for key, value in gap_dic.items()]
    return gap_list

def length_of_lis(nums):
    if not nums:
        return 0
    d = []
    for num in nums:
        pos = bisect_left(d, num)
        if pos == len(d):
            d.append(num)
        else:
            d[pos] = num
    return len(d)/len(nums)

def longest_common_subsequence(list1, list2):
    list1_dict = {key: value for key, value in list1.items()}
    dp = [list1_dict[element] for element in list2 if element in list1_dict]
    return length_of_lis(dp)

def create_position_dict(data):
    position_dict = defaultdict(lambda: 0)
    for content, positions in data.items():
        for pos_tuple in positions:
            if pos_tuple:
                position_dict[pos_tuple[0]] = content

    return position_dict

def create_boundaries(intervals):
    boundaries = [0]
    for start, end in intervals:
        if not boundaries or start > boundaries[-1]:
            boundaries.append(start)
        if end > boundaries[-1]:
            boundaries.append(end)
    boundaries.append(1000000000)
    return [(boundaries[i], boundaries[i + 1]) for i in range(0, len(boundaries) - 1, 2)]

def assign_intervals(new_intervals, original_intervals):
    boundaries = create_boundaries(original_intervals)
    result=[]
    start=new_intervals[0][0]
    end=new_intervals[0][1]
    assigned = False
    for i, (boundary_start, boundary_end) in enumerate(boundaries):
        if boundary_start <= start < boundary_end and boundary_start < end <= boundary_end:
            result.append(i)
            assigned = True
            break
    if not assigned:
        result.append(None)
    return result[0]

def integrateReadLine(fa, prefix): 
    file_name = fa.split("/")[-1]
    fi = open(fa, "r")
    lines = fi.readlines()
    fi.close()
    fo = open(prefix + "_" + file_name, "w")
    n = "1"
    for line in lines:
        if line[0] == ">":
            if n == "1":
                fo.write(line)
            else:
                fo.write("\n" + line)
        else:
            n = "larger than 1"
            line = line.strip()
            fo.write(line)
    fo.close()
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) + ", " + "Integrate read lines in " + fa + " file done.")

def reverseCompleKmer(kmer):
    reverse_complementary_kmer_list = []
    base_dic = {"A":"T","T":"A","C":"G","G":"C","N":"N","a":"T","t":"A","c":"G","g":"C","n":"N"}
    for base in list(kmer[::-1]):
        reverse_complementary_kmer_list.append(base_dic[base])
    reverse_complementary_kmer = "".join(reverse_complementary_kmer_list)
    return min([kmer.upper(), reverse_complementary_kmer.upper()])

def refUniqueKmerSearch(ref, kmer_length):
    fi = open(ref, "r")
    lines = fi.readlines()
    fi.close()
    unique_kmer = {}
    sumKmer_dic = {}
    for line in lines:
        line = line.strip()
        if line[0] != ">":
            for i in range(len(line) - kmer_length + 1):
                reverse_comple_kmer = reverseCompleKmer(line[i:i + kmer_length])
                if reverse_comple_kmer in sumKmer_dic:
                    x=unique_kmer.pop(reverse_comple_kmer, None)
                else:
                    unique_kmer[reverse_comple_kmer] = 0
                    sumKmer_dic[reverse_comple_kmer] = 0
    sumKmer_dic =None
    internal_refChr_uniKmer = {}
    refChr_uniKmer_Pos = {}
    chr_name = ""
    for line in lines:
        line = line.strip().split()[0]
        if line != "":
            if line[0] == ">":
                chr_name = line[1:]
                internal_refChr_uniKmer = {}
            else:
                for i in range(len(line) - kmer_length + 1):
                    kmer = line[i:i + kmer_length]
                    kmer = reverseCompleKmer(kmer)
                    if kmer in unique_kmer.keys():

                        internal_refChr_uniKmer[(i + 1, i + kmer_length)]=0

                refChr_uniKmer_Pos[chr_name] = internal_refChr_uniKmer

    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) + ", " + "Unique kmer count in " + ref + " file done.")
    return refChr_uniKmer_Pos

def getRefUniKmerPos(ref, kmer_length):
    eachchr_unikmer = refUniqueKmerSearch(ref, kmer_length)[0]
    data = open(ref, "r")
    lines = data.readlines()
    data.close()
    kmerpos_dic = {}
    chr_name = ""
    for line in lines:
        line = line.strip().split()[0]
        if line[0] == ">":
            chr_name = line[1:]
            kmerpos_dic[chr_name] = {}
    for line in lines:
        line = line.strip().split()[0]
        if line[0] == ">":
            chr_name = line[1:]
        else:
            for i in range(len(line) - kmer_length + 1):
                if reverseCompleKmer(line[i:i + kmer_length]) in eachchr_unikmer[chr_name].keys():
                    kmerpos_dic[chr_name][(i + 1, i + kmer_length)] = 0
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) + ", " + "Get reference unique kmer position done.")
    return kmerpos_dic

def removeOverlapKmer(classed_kmer_dic):
    set_option('display.max_columns', None)
    set_option('display.max_rows', None)
    external_dic = {}
    internal_dic = {}
    for key, value in classed_kmer_dic.items():
        internal_dic = {}
        value_list = list(value.items())
        value_list.sort(key=lambda x: x[0][0], reverse=False)
        mydic = {}
        mydic["start_pos"] = [value_list[0][0][0]]
        mydic["end_pos"] = [value_list[0][0][1]]
        mydic["front_end_pos"] = [0]
        last_end_pos = value_list[0][0][1]

        for i in range(1,len(value_list)):
            mydic["start_pos"].append(value_list[i][0][0])
            mydic["end_pos"].append(value_list[i][0][1])
            mydic["front_end_pos"].append(last_end_pos)
            last_end_pos = value_list[i][0][1]
        mydataframe = DataFrame(mydic)
        mydataframe = mydataframe[(mydataframe["front_end_pos"] - mydataframe["start_pos"]) < 0]
        for i in mydataframe.index:
            internal_dic[(mydataframe["start_pos"][i],mydataframe["end_pos"][i])] = 0
        external_dic[key] = internal_dic
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ", " + "Remove overlaps between kmers done.")
    return external_dic

def getHeaderKmer(fasta, headerkmerpos_dic, kmer_length):
    data = open(fasta, "r")
    lines = data.readlines()
    data.close()
    kmer_pos_dic = {}
    seqname_kmer_dic = {}
    sequence_name = ""
    for line in lines:
        line = line.strip().split()[0]
        if line[0] == ">":
            sequence_name = line[1:]
            kmer_pos_dic = {}
        else:
            for key, value in headerkmerpos_dic.items():
                if sequence_name == key:
                    for pos in value.keys():
                        kmer_pos_dic[reverseCompleKmer(line[pos[0] - 1:pos[0] + kmer_length - 1])] = [pos]
                    seqname_kmer_dic[sequence_name] = kmer_pos_dic
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ", " + "Get header kmer done.")
    return seqname_kmer_dic

def randomselectkmer(refkmerdic,sample_num):
    kmerdic = {}
    for key,value in refkmerdic.items():
        for k,v in value.items():
            kmerdic[k] = 0
    seed(111)
    mylist = sample(list(kmerdic.items()),sample_num)
    newkmerdic = {}
    for i in mylist:
        newkmerdic[i[0]] = 0
    external_dic = {}
    for key,value in refkmerdic.items():
        internal_dic = {}
        for k,v in value.items():
            if k in newkmerdic.keys():
                internal_dic[k] = v
        external_dic[key] = internal_dic
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ", random select 200000 reference kmer done.")
    return external_dic

def getAsbkmer(scafile,reffile,refkmer_dic,kmer_length):
    refkmer = {}
    for key,value in refkmer_dic.items():
        for k,v in value.items():
            refkmer[k] = 0
    file = open(scafile,"r")
    lines = file.readlines()
    file.close()
    external_dic = {}
    internal_dic = {}
    scaname = ""
    for line in lines:
        line = line.strip().split()[0]
        if line[0] == ">":
            scaname = line[1:]
            internal_dic = {}
        else:
            for i in range(len(line) - kmer_length + 1):
                if reverseCompleKmer(line[i:i + kmer_length]) in refkmer.keys():
                    internal_dic[reverseCompleKmer(line[i:i + kmer_length])] = internal_dic.get(reverseCompleKmer(line[i:i + kmer_length]),[]) + [(i + 1, i + kmer_length)]
            external_dic[scaname] = internal_dic
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) + ", " + "Unique kmer from " + reffile + " file count in " + scafile + " file done.")
    return external_dic

def assembleAssessment(ref_fa, asb_fa, kmer_length, sample_num, prefix):
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ", Begin.")
    integrateReadLine(ref_fa, prefix)
    integrateReadLine(asb_fa, prefix)
    refUniKmerPos = refUniqueKmerSearch(prefix + "_" + ref_fa.split("/")[-1], kmer_length)
    refHeaderKmerPos = removeOverlapKmer(refUniKmerPos)
    refUniKmerPos = None
    if sample_num == "all":
        refHeaderKmer = getHeaderKmer(prefix + "_" + ref_fa.split("/")[-1], refHeaderKmerPos, kmer_length)
    else:
        refheaderKmer = getHeaderKmer(prefix + "_" + ref_fa.split("/")[-1], refHeaderKmerPos, kmer_length)
        refHeaderKmer = randomselectkmer(refheaderKmer,int(sample_num))
    refheaderkmer = None
    refHeaderKmerPos = None
    asbHeaderKmer = getAsbkmer(prefix + "_" + asb_fa.split("/")[-1],prefix + "_" + ref_fa.split("/")[-1],refHeaderKmer,kmer_length)
    with  open("random_refHeaderkmer.dic","w") as file:

        file.write(str(refHeaderKmer))
    file.close()
    with open("refHeaderkmer.out", "w") as file:
        file.write("Chr\tKmer\tStart\tEnd\n")
        for key, value in refHeaderKmer.items():
            for k, v in value.items():
                file.write("{}\t{}\t{}\t{}\n".format(key, k, v[0][0], v[0][1]))
            file.write("\n")
    file.close()
    with open("asbHeaderkmer.out", "w") as file:
        file.write("scaffold\tKmer\tStart\tEnd\n")
        for key, value in asbHeaderKmer.items():
            for k, v in value.items():
                for i in v:
                    file.write("{}\t{}\t{}\t{}\n".format(key, k, i[0], i[1]))
            file.write("\n")
    file.close()
    asbKmerCount = {}
    for Key, Value in asbHeaderKmer.items():
        for k, v in Value.items():
            asbKmerCount[k] = asbKmerCount.get(k, 0) + len(v)
    singleCopyKmerNum = 0
    duplicateKmerNum = 0
    colKmerNum = 0
    for key, value in refHeaderKmer.items():
        for k, v in value.items():
            colKmerNum += 1
    for key, value in asbKmerCount.items():
        if value == 1:
            singleCopyKmerNum += 1
        else:
            duplicateKmerNum += 1
    singleCopy = singleCopyKmerNum / colKmerNum
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ", " + "Compute single copy rate done.")
    duplicateRate = duplicateKmerNum / colKmerNum
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()) + ", " + "Compute duplicate rate done.")
    scakmer_dic = {}
    for key, value in asbHeaderKmer.items():
        scakmer_dic[key] = {}
    for key, value in asbHeaderKmer.items():
        for k, v in value.items():
            if k not in scakmer_dic[key].keys():
                scakmer_dic[key][k] = 0
    chrkmer_dic = {}
    for key, value in refHeaderKmer.items():
        chrkmer_dic[key] = {}
    for key, value in refHeaderKmer.items():
        for k, v in value.items():
            if k not in chrkmer_dic[key].keys():
                chrkmer_dic[key][k] = 0
    chr_to_kmerlist = {}
    scaffold_to_chr = {}
    for sca, scakmer in asbHeaderKmer.items():
        chr_to_kmerlist = {}
        for chr, chrkmer in refHeaderKmer.items():
            chr_to_kmerlist[chr] = {}
        scaffold_to_chr[sca] = chr_to_kmerlist
    for asb_key, asb_value in scakmer_dic.items():
        for i in asb_value.keys():
            for ref_key, ref_value in chrkmer_dic.items():
                if i in ref_value.keys() and i not in scaffold_to_chr[asb_key][ref_key].keys():
                    scaffold_to_chr[asb_key][ref_key][i] = 0
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) + ", " + "Get collective header kmer from reference and assemble file done.")
    scaffold_to_kmer = {}
    for key, value in scaffold_to_chr.items():
        kmerlist = {}
        for k, v in value.items():
            for kmer, it in v.items():
                kmerlist[kmer] = 0
        scaffold_to_kmer[key] = kmerlist
    scakmer_dic = None
    chrkmer_dic = None
    proportion_of_the_largest_categories = 0
    largest_categories_ex = {}
    largest_categories_in = {}
    numerator = 0
    denominator = 0
    for key, value in scaffold_to_chr.items():  
        largest_categories_num = 0
        largest_categories_in = {}
        for k, v in value.items():
            if len(v) > largest_categories_num:
                largest_categories_num = len(v)
        for k, v in value.items():
            if len(v) == largest_categories_num:
                largest_categories_in[k] = v
                largest_categories_ex[key] = largest_categories_in
        numerator += largest_categories_num
        denominator += len(scaffold_to_kmer[key])
    proportion_of_the_largest_categories += float(numerator) / float(denominator)
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) + ", " + "Compute proportion of the largest categories done.")
    aveEachScaDistance_sum = 0
    scaffold_to_chr = None
    with open("largest_categories.dic","w") as fi:
        fi.write(str(largest_categories_ex))
    fi.close()
    with open("refKmer.dic","w") as fi:
        fi.write(str(refHeaderKmer))
    fi.close()
    with open("asbKmer.dic","w") as fi:
        fi.write(str(asbHeaderKmer))
    fi.close()
    scaffold_nums = 0
    LIS = 0
    ti = 0
    tj = 0
    aveeacsacffolding = 0
    aveeacsacffolding1 = 0
    gap_list=find_gaps(prefix + "_" + asb_fa.split("/")[-1])
    for key, value in largest_categories_ex.items():
        ref_base_pos_dic = {}
        asb_base_pos_dic = {}
        asb_ref_part_dic = {}
        t = 0
        eachScaDistance_sum = 0
        for k, v in value.items():
            if len(v) > 1:
                scaffold_nums += 1
                for i in range(len(v)):
                    asb_ref_base_pos_dic = {}
                    ref_base_pos_dic[list(v.items())[i][0]] = refHeaderKmer[k][list(v.items())[i][0]]
                    asb_base_pos_dic[list(v.items())[i][0]] = asbHeaderKmer[key][list(v.items())[i][0]]
                    for xi, xj in gap_list:
                        if xi == key and xj != []:
                            jud_factor = assign_intervals(asbHeaderKmer[key][list(v.items())[i][0]], xj)
                            if jud_factor != None:
                                asb_ref_base_pos_dic[list(v.items())[i][0]] = asbHeaderKmer[key][list(v.items())[i][0]], \
                                                                              refHeaderKmer[k][list(v.items())[i][0]]
                                if jud_factor in asb_ref_part_dic:
                                    asb_ref_part_dic[jud_factor].update(asb_ref_base_pos_dic)
                                else:
                                    asb_ref_part_dic[jud_factor] = asb_ref_base_pos_dic
                ref_base_list = list(ref_base_pos_dic.items())
                ref_base_list.sort(key=lambda x: x[1][0][0], reverse=False)
                for i in range(len(ref_base_list) - 1):
                    refKmerDistance = abs(ref_base_list[i + 1][1][0][0] - ref_base_list[i][1][0][0])
                    asbKmerDistance = {}
                    for j in range(len(asb_base_pos_dic[ref_base_list[i][0]])):
                        for m in range(len(asb_base_pos_dic[ref_base_list[i + 1][0]])):
                            distance = abs(asb_base_pos_dic[ref_base_list[i][0]][j][0] - asb_base_pos_dic[ref_base_list[i + 1][0]][m][0])
                            asbKmerDistance[distance] = 0
                    for j in asbKmerDistance.keys():
                        eachScaDistance_sum += abs(j - refKmerDistance)
                        t += 1
                aveEachScaDistance_sum += eachScaDistance_sum / t
                ref_notreverse_list = list(ref_base_pos_dic.items())
                ref_reverse_list = list(ref_base_pos_dic.items())
                ref_notreverse_list.sort(key=lambda x: x[1][0][0], reverse=False)
                ref_reverse_list.sort(key=lambda x: x[1][0][0], reverse=True)
                position_dict = create_position_dict(asb_base_pos_dic)
                sorted_items = sorted(position_dict.items())
                sorted_values = [value for key, value in sorted_items]
                result = {item[0]: index + 1 for index, item in enumerate(ref_notreverse_list)}
                result_revaerse = {item[0]: index + 1 for index, item in enumerate(ref_reverse_list)}
                LIS += max(longest_common_subsequence(result, sorted_values),
                           longest_common_subsequence(result_revaerse, sorted_values))
        eachsac_distance = 0
        count = 0
        for pos, area_dic in asb_ref_part_dic.items():
            dif_sum = 0
            if pos + 1 in asb_ref_part_dic:

                for f in range(min(len(asb_ref_part_dic[pos].items()), len(asb_ref_part_dic[pos + 1].items()))):
                    area_i_list = ([(key, value) for key, value in asb_ref_part_dic[pos].items()])
                    area_iadd1_list = ([(key, value) for key, value in asb_ref_part_dic[pos + 1].items()])
                    dif_sum += abs(abs((area_iadd1_list[f][1][0][0][0] - area_i_list[f][1][0][0][0])) - abs(
                        (area_iadd1_list[f][1][1][0][0] - area_i_list[f][1][1][0][0])))
                count += 1
            else:
                break
            avg_sum = dif_sum / min(len(asb_ref_part_dic[pos].items()), len(asb_ref_part_dic[pos + 1].items()))
            eachsac_distance += avg_sum
        if count != 0:
            ti=ti+1
        eachsac_distance = eachsac_distance / count if count else 0
        aveeacsacffolding = aveeacsacffolding + eachsac_distance
        eachsac_distance1 = 0
        count1 = 0
        keys_list= [key for key in asb_ref_part_dic]
        for i in range(len(keys_list)):

            try:
                keys_list[i]
            except IndexError :
                pos=-1
            else :
                pos=keys_list[i]
            try:
                keys_list[i+1]
            except IndexError :
                pos1=-1
            else :
                pos1=keys_list[i+1]
            dif_sum1 = 0
            if pos1 in asb_ref_part_dic:
                for f in range(min(len(asb_ref_part_dic[pos].items()), len(asb_ref_part_dic[pos1].items()))):
                    area_i_list = ([(key, value) for key, value in asb_ref_part_dic[pos].items()])
                    area_iadd1_list = ([(key, value) for key, value in asb_ref_part_dic[pos1].items()])
                    dif_sum1 += abs(abs((area_iadd1_list[f][1][0][0][0] - area_i_list[f][1][0][0][0])) - abs(
                        (area_iadd1_list[f][1][1][0][0] - area_i_list[f][1][1][0][0])))
                count1 += 1
            else:
                break
            avg_sum1 = dif_sum1 / min(len(asb_ref_part_dic[pos].items()), len(asb_ref_part_dic[pos1].items()))
            eachsac_distance1 += avg_sum1
        if count1 != 0:
            tj=tj+1
        eachsac_distance1 = eachsac_distance1 / count1 if count1 else 0
        aveeacsacffolding1 = aveeacsacffolding1 + eachsac_distance1
    ave_distance_diff = aveEachScaDistance_sum / scaffold_nums
    LIS_distance_diff = LIS / scaffold_nums
    if ti !=0:
        ave_scaffold_diff = aveeacsacffolding / ti
    else:
        ave_scaffold_diff = aveeacsacffolding1 / tj if tj else 0
    print("[INFO] " + time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()) + ", " + "Compute average distance difference done." + "\n")
    print("{:<24.7f}Complete\n"
          "{:<24.7f}Complete and single-copy\n"
          "{:<24.7f}Complete and duplicated\n"
          "{:<24.7f}Proportion of the largest categories\n"
          "{:<24.7f}ave distance diff\n"
          "{:<24.7f}Inter-contig Average Distance Difference\n"
          "{:<24.7f}Length of the k-mer-based Longest Increasing Subsequence\n".format(singleCopy + duplicateRate, singleCopy, duplicateRate,proportion_of_the_largest_categories, ave_distance_diff,LIS_distance_diff,ave_scaffold_diff))
    with open("result.report", "w") as fi:
        fi.write("{:<24.7f}Complete\n"
                 "{:<24.7f}Complete and single-copy\n"
                 "{:<24.7f}Complete and duplicated\n"
                 "{:<24.7f}Proportion of the largest categories\n"
                 "{:<24.7f}ave distance diff\n"
                 "{:<24.7f}Inter-contig Average Distance Difference\n"
                 "{:<24.7f}Length of the k-mer-based Longest Increasing Subsequence\n".format(singleCopy + duplicateRate, singleCopy, duplicateRate,proportion_of_the_largest_categories, ave_distance_diff,LIS_distance_diff,ave_scaffold_diff))
        fi.close()

if __name__ == "__main__":
    opt = ArgParse()
    path = os.getcwd()
    asb_seq = opt.input
    ref_seq = opt.reference
    kmer_len = opt.kmer_length
    prefix = opt.out_prefix
    if opt.sample != "all":
        sample_num = int(opt.sample)
    else:
        sample_num = opt.sample
    asb_seq_path = os.path.join(path,asb_seq)
    ref_seq_path = os.path.join(path,ref_seq)
    assembleAssessment(ref_seq_path, asb_seq_path, kmer_len, sample_num, prefix)
