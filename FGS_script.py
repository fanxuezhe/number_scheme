from hmmer_parser import *
#2018.12.10  FGS_多样性分析
#正规一代测序运行流程scFv
#第一步利用FGS模块将核苷酸序列合并。并且拆分为VH VL 两个文件，同时记录下来含有终止密码子的序列
#提前定义好特异性的参数

F=["ATTTGGCACAAGCGGCCGCT","ATAAACGCGAGGCCGAGGCG","AAAGAGAGGCCGAGTCGGCC","CAGTTTTGGCACAGGCGGCC","CAGTTTTGGCACGGGCGGCC","AAAGAGAGGCCGAGGCGGCC","TTTTGGCACAAGCGGCCGCT","CAGTTTTGGCACACGCGGCC","TAGTTTTGGCACAGGCGGCC","AACGCGAGGCCGAGGCGGCC","CAGTTTTGGCCCAGGCGGCC","TTTTGGCACAAGCGGCCGCT","TTTTGGCCCAAGCGGCCGCT","CAGAATTGGCCCAGGCGGCC"]
R=["TCACAAGGGCCATCGTGTTC","ACTGTGGCGGCCCCATCTGT","GGCCGTCGTGGTCGACCTGC","GGGCCGTCGTGGTCGACCTG","GGGCCGTTGGTGGTCGACCT","GGGCCGTCGGTGGTCGACCT","GGGCCGTCGGTGGTCGACTG","GGGCCGTCGGTGGCGACTGC","GGGCCGTCGGTGGTCAACTG","AGGCCGTCGGTGGTCGACTG","GGCCGTCGGTGGTCGACCTG","GGCCCGTCGGTGGTCGACCT","GCCTCCACAAGGGCCATCGG","CAGTTTTGGCCCAGGCGGCC","GCCTCCACCAAGGGGCCGTC","GGGCCGTCGGTGTTCCCCCT","GCCTCCACCAAGGGGCCGTC","GCCTCCACCAAGGGCCCATC","GCCTCCACCAAGGGGCCATC"]
split_all=["RAGGGGSGGGGSGGLGGAAV","RAGGGGSGGGGSGGLGCAAA","RAGGGGSGGGGSGGLGGAGA","RAGGGGSGGGGSGGLCGAAA","RAGGGGSGGGGLGGLGGAAA","RAGGGGSGGGESGGLGGAAA","RAGGGGSGGGSGGLGGAAA","RAGGGGSGGGGSGGLGGTAA","APGGGGSGGGGSGGLGGAAA","RAEGRGSGGGGSGGLGGAAA","RAGEGGSGGGGSGGLGGAAA","APEEEGSGGGGSGGLGGAAA","RAGGEGSGGGGGGAAA","RAGGGDSGGGGSGGLGGAAA","RAGGGGSGGGGSGGLGGAAA","APGGGGSGGGGSGGLGGAA","APGGGGSGGGGSGRTGGAAA","RAGGGGSGGGGGGAAA","RAGGGSGGGGSGGLGGAAA","RAGGGCSGGGGSGGLGGAAA","RAGGEGSGGGGSGGLGGAAA","RAVGGGSGGGGSGGLGGAAA","RAGRGGSGGGGSGELGGAAA","RAGGGGSGGSGSGGLGGAAA","RAGGGSSGGGGSGGLGGAAA","RAGGGGSEGGGSGGLGGAAA"]#这里需要在循环过程中不断地添加新的linker序列
split_all_short=["RAGGGGSGGG"]
#正式进行分析
#首先新建文件夹来储存所有的轻重链的fasta文件
import shutil
#首先定义将用来储存fasta文件的文件夹所在位置以及异常序列文件夹所在的位置，这一部分运行完后及时注释掉，以免影响后续的数据分析
file_repo=r"E:\Program\project\yanli\FGS_2018.12.11\all_files"#这个文件夹是所有含有测序数据的文件夹所在的文件夹
#定义每一个含有一代测序文件的文件夹，并且对其进行分析
#首先定义一个总的目录，在这里所有含有一代测序文件的文件夹都包含在其中
result_folder=(file_repo+"\\result")
all_fasta=(result_folder+"\\fasta_file")
all_folders=get_folder(file_repo)
print all_folders
#organize_result(file_repo,all_folders)#该步骤运行完之后及时注释，防止覆盖原来的文件
#scfv_analysis(file_repo,all_folders)

#########################第二步###########################################
#接下来进行HMMSCAN分析，并且将所有的比对结果和fasta文件都拷贝到一个文件夹中
###########################################################################



#############################第三步###########################################
##############################################################################
######接下来是不用修改的部分（将所有需要新建的文件夹在这里建立，之后假如有需要修改的话相应的软件中也需要修改）########################
match_repo=(file_repo+"\\all_match_file")#这一行不用注释
VH_folder=(match_repo+"\\VH")#这一行不用注释
VL_folder=(match_repo+"\\VL")#这一行不用注释
hmm_out_folder=(file_repo+"\\all_hmm_out")#这个文件夹需要从其他地方拷贝过来然后解压在file_repo文件夹中，该行不用注释
VH_VL_folder=(match_repo+"\\VH_VL_file")#这一行不用注释
VH_CDR3=(VH_folder+"\\VH_CDR3")#这个文件夹将会进行CDR3的频率分析
C_file=(match_repo+"\\C_file")#这个文件夹将会储存所有CDR以及CDR3的分析
C_combined=(C_file+"\\VH_VL_file")#这个文件夹将会用来分析C文件的合并频率
combined_folder=(VH_VL_folder+"\\combined_VH_VL")#这一行不用注释
rm_folder=(C_combined+"\\rm_folder")#这一行不用注释
C_combined_all=(C_combined+"\\all_combined")
combined_all=(combined_folder+"\\all_combined")
C_CDR3=(C_file+"\\C_CDR3")
rm_folder_CDR3=(C_CDR3+"\\without_species")
species_CDR3_only_folder=(C_CDR3+"\\CDR3_only")
CDR3_only_folder=(VH_CDR3+"\\CDR3_only")
species_CDR3_formal_id_folder=(C_CDR3+"\\formal_id")
CDR3_formal_id_folder=(VH_CDR3+"\\formal_id")
# 

mkdir_folder(match_repo)
mkdir_folder(VH_folder)
mkdir_folder(VL_folder)
mkdir_folder(VH_VL_folder)
mkdir_folder(VH_CDR3)
mkdir_folder(C_file)
mkdir_folder(C_combined)
mkdir_folder(combined_folder)
mkdir_folder(rm_folder)
mkdir_folder(C_combined_all)
mkdir_folder(combined_all)
mkdir_folder(C_CDR3)
mkdir_folder(rm_folder_CDR3)
mkdir_folder(species_CDR3_only_folder)
mkdir_folder(CDR3_only_folder)#这里运行一遍以后及时注释，否则会多次对文件进行提取CDR3的操作
mkdir_folder(species_CDR3_formal_id_folder)
mkdir_folder(CDR3_formal_id_folder)



#1.取出CDR序列，随后将所有的hmm输出文件以及fasta文件拷贝到match_repo文件夹中，并且解压
#然后将所有的轻链和重链的序列和hmmscan输出文件分别拷贝到VH和VL文件夹中
type_char="VH"
mkdir="F"
cp_all_file(all_fasta,VH_folder,type_char,mkdir)
cp_all_file(hmm_out_folder,VH_folder,type_char,mkdir)
type_char="VL"
mkdir="F"
cp_all_file(all_fasta,VL_folder,type_char,mkdir)
cp_all_file(hmm_out_folder,VL_folder,type_char,mkdir)
#2.分别对VH和VL进行分析，首先分析VH
scheme="Longest_CDR_H"
dict_match=match_fasta_txt(VH_folder)
print dict_match
get_cdr_from_hmmer(VH_folder,dict_match,scheme)
scheme="Longest_CDR_L"
dict_match=match_fasta_txt(VL_folder)
print dict_match
get_cdr_from_hmmer(VL_folder,dict_match,scheme)
#3.取出CDR序列文件，分别放置在VH_VL以及VH单独的文件夹中，便于后续的使用    
type_char="CDR_sequences"
mkdir="F"
cp_all_file(VH_folder,VH_VL_folder,type_char,mkdir)
cp_all_file(VL_folder,VH_VL_folder,type_char,mkdir)
cp_all_file(VH_folder,VH_CDR3,type_char,mkdir)
#################################################第四步########################################################
#######利用merge_VH_VL将VH和VL的序列的CDR区域合并
#1.首先将重链和轻链的序列拷贝到combined_folder中去，便于后续的合并步骤的进行
pos_char=".fasta"#将全部的fasta文件拷贝到新的文件夹中去
pos="end"
cp_start_end(VH_VL_folder,combined_folder,pos_char,pos)
pos_char="C"#这里将C开头的文件拷贝到未来需要去除物种的文件夹中去
pos="start"
cp_start_end(VH_VL_folder,C_combined,pos_char,pos)

#2.然后进行合并（VH+VL）先处理所有的数据，在处理C开头的数据
#首先将C开头的文件进行处理，后续的操作和常规的一样，便于程序的合并
species="rabbit"
divide_species_cdr(C_combined,species)
type_char=("without")
mkdir="F"
cp_all_file(C_combined,rm_folder,type_char,mkdir)
#对文件的ID进行处理以保证在合并的时候不会出现问题
rm_file="T"
split_char="_"
pos_start=0
pos_end=-3
split_id(combined_folder,split_char,pos_start,pos_end,rm_file)
merge_VH_VL(combined_folder)
split_id(rm_folder,split_char,pos_start,pos_end,rm_file)
merge_VH_VL(rm_folder)
#新建文件夹储存C开头的VH_VL合并文件以及所有合并文件
type_char="combined_VH_VL_CDR"
mkdir="F"
cp_all_file(rm_folder,C_combined_all,type_char,mkdir)
cp_all_file(combined_folder,combined_all,type_char,mkdir)

############################################第五步#########################################
###############################################################################################
#########################将所有重链的CDR3的序列提取出来#################################
#1.建立储存C开头的CDR3的文件夹，并且拷贝进去C开头的CDR3

#2.将C开头的重链序列转移到C_CDR3中
type_char="C"
mkdir="F"
cp_all_file(VH_CDR3,C_CDR3,type_char,mkdir)
#3.按照物种进行分析
species="rabbit"
divide_species_cdr(C_CDR3,species)
#4.新建文件夹储存不含特定物种的重链CDR3
type_char="without"
mkdir="F"
cp_all_file(C_CDR3,rm_folder_CDR3,type_char,mkdir)
#获取CDR3,并且储存在CDR3_only文件夹中（包括总的CDR3和只是C开头的）
get_cdr3(rm_folder_CDR3)
get_cdr3(VH_CDR3)
type_char="CDR3_only"
mkdir="F"
cp_all_file(rm_folder_CDR3,species_CDR3_only_folder,type_char,mkdir)
cp_all_file(VH_CDR3,CDR3_only_folder,type_char,mkdir)
#编辑ID
rm_file="T"
split_char="_"
pos_start=0
pos_end=-3
split_id(species_CDR3_only_folder,split_char,pos_start,pos_end,rm_file)
split_id(CDR3_only_folder,split_char,pos_start,pos_end,rm_file)
#拷贝编辑ID后的序列拷贝到新的文件夹中
type_char="edited_id"
mkdir="F"
cp_all_file(species_CDR3_only_folder,species_CDR3_formal_id_folder,type_char,mkdir)
cp_all_file(CDR3_only_folder,CDR3_formal_id_folder,type_char,mkdir)
#对所有的需要进行频率分析的文件夹进行频率分析
#第一部分是对常规序列的分析
all_path=[species_CDR3_formal_id_folder,CDR3_formal_id_folder,C_combined_all,combined_all]
for path_n in all_path:
    if path_n.find("CDR3")!=-1:
        len_out=1
        freq_analysis_fasta_ID(path_n,len_out)
    elif path_n.find("all_combined")!=-1:
        len_out=6
        freq_analysis_fasta_ID(path_n,len_out)
#将对应的文件拷贝到相应的文件夹中
cp_file_to_result(file_repo)
