from hmmer_parser import *
#2018.12.10  FGS_�����Է���
#����һ��������������scFv
#��һ������FGSģ�齫���������кϲ������Ҳ��ΪVH VL �����ļ���ͬʱ��¼����������ֹ�����ӵ�����
#��ǰ����������ԵĲ���

F=["ATTTGGCACAAGCGGCCGCT","ATAAACGCGAGGCCGAGGCG","AAAGAGAGGCCGAGTCGGCC","CAGTTTTGGCACAGGCGGCC","CAGTTTTGGCACGGGCGGCC","AAAGAGAGGCCGAGGCGGCC","TTTTGGCACAAGCGGCCGCT","CAGTTTTGGCACACGCGGCC","TAGTTTTGGCACAGGCGGCC","AACGCGAGGCCGAGGCGGCC","CAGTTTTGGCCCAGGCGGCC","TTTTGGCACAAGCGGCCGCT","TTTTGGCCCAAGCGGCCGCT","CAGAATTGGCCCAGGCGGCC"]
R=["TCACAAGGGCCATCGTGTTC","ACTGTGGCGGCCCCATCTGT","GGCCGTCGTGGTCGACCTGC","GGGCCGTCGTGGTCGACCTG","GGGCCGTTGGTGGTCGACCT","GGGCCGTCGGTGGTCGACCT","GGGCCGTCGGTGGTCGACTG","GGGCCGTCGGTGGCGACTGC","GGGCCGTCGGTGGTCAACTG","AGGCCGTCGGTGGTCGACTG","GGCCGTCGGTGGTCGACCTG","GGCCCGTCGGTGGTCGACCT","GCCTCCACAAGGGCCATCGG","CAGTTTTGGCCCAGGCGGCC","GCCTCCACCAAGGGGCCGTC","GGGCCGTCGGTGTTCCCCCT","GCCTCCACCAAGGGGCCGTC","GCCTCCACCAAGGGCCCATC","GCCTCCACCAAGGGGCCATC"]
split_all=["RAGGGGSGGGGSGGLGGAAV","RAGGGGSGGGGSGGLGCAAA","RAGGGGSGGGGSGGLGGAGA","RAGGGGSGGGGSGGLCGAAA","RAGGGGSGGGGLGGLGGAAA","RAGGGGSGGGESGGLGGAAA","RAGGGGSGGGSGGLGGAAA","RAGGGGSGGGGSGGLGGTAA","APGGGGSGGGGSGGLGGAAA","RAEGRGSGGGGSGGLGGAAA","RAGEGGSGGGGSGGLGGAAA","APEEEGSGGGGSGGLGGAAA","RAGGEGSGGGGGGAAA","RAGGGDSGGGGSGGLGGAAA","RAGGGGSGGGGSGGLGGAAA","APGGGGSGGGGSGGLGGAA","APGGGGSGGGGSGRTGGAAA","RAGGGGSGGGGGGAAA","RAGGGSGGGGSGGLGGAAA","RAGGGCSGGGGSGGLGGAAA","RAGGEGSGGGGSGGLGGAAA","RAVGGGSGGGGSGGLGGAAA","RAGRGGSGGGGSGELGGAAA","RAGGGGSGGSGSGGLGGAAA","RAGGGSSGGGGSGGLGGAAA","RAGGGGSEGGGSGGLGGAAA"]#������Ҫ��ѭ�������в��ϵ�����µ�linker����
split_all_short=["RAGGGGSGGG"]
#��ʽ���з���
#�����½��ļ������������е���������fasta�ļ�
import shutil
#���ȶ��彫��������fasta�ļ����ļ�������λ���Լ��쳣�����ļ������ڵ�λ�ã���һ�����������ʱע�͵�������Ӱ����������ݷ���
file_repo=r"E:\Program\project\yanli\FGS_2018.12.11\all_files"#����ļ��������к��в������ݵ��ļ������ڵ��ļ���
#����ÿһ������һ�������ļ����ļ��У����Ҷ�����з���
#���ȶ���һ���ܵ�Ŀ¼�����������к���һ�������ļ����ļ��ж�����������
result_folder=(file_repo+"\\result")
all_fasta=(result_folder+"\\fasta_file")
all_folders=get_folder(file_repo)
print all_folders
#organize_result(file_repo,all_folders)#�ò���������֮��ʱע�ͣ���ֹ����ԭ�����ļ�
#scfv_analysis(file_repo,all_folders)

#########################�ڶ���###########################################
#����������HMMSCAN���������ҽ����еıȶԽ����fasta�ļ���������һ���ļ�����
###########################################################################



#############################������###########################################
##############################################################################
######�������ǲ����޸ĵĲ��֣���������Ҫ�½����ļ��������ｨ����֮���������Ҫ�޸ĵĻ���Ӧ�������Ҳ��Ҫ�޸ģ�########################
match_repo=(file_repo+"\\all_match_file")#��һ�в���ע��
VH_folder=(match_repo+"\\VH")#��һ�в���ע��
VL_folder=(match_repo+"\\VL")#��һ�в���ע��
hmm_out_folder=(file_repo+"\\all_hmm_out")#����ļ�����Ҫ�������ط���������Ȼ���ѹ��file_repo�ļ����У����в���ע��
VH_VL_folder=(match_repo+"\\VH_VL_file")#��һ�в���ע��
VH_CDR3=(VH_folder+"\\VH_CDR3")#����ļ��н������CDR3��Ƶ�ʷ���
C_file=(match_repo+"\\C_file")#����ļ��н��ᴢ������CDR�Լ�CDR3�ķ���
C_combined=(C_file+"\\VH_VL_file")#����ļ��н�����������C�ļ��ĺϲ�Ƶ��
combined_folder=(VH_VL_folder+"\\combined_VH_VL")#��һ�в���ע��
rm_folder=(C_combined+"\\rm_folder")#��һ�в���ע��
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
mkdir_folder(CDR3_only_folder)#��������һ���Ժ�ʱע�ͣ�������ζ��ļ�������ȡCDR3�Ĳ���
mkdir_folder(species_CDR3_formal_id_folder)
mkdir_folder(CDR3_formal_id_folder)



#1.ȡ��CDR���У�������е�hmm����ļ��Լ�fasta�ļ�������match_repo�ļ����У����ҽ�ѹ
#Ȼ�����е����������������к�hmmscan����ļ��ֱ𿽱���VH��VL�ļ�����
type_char="VH"
mkdir="F"
cp_all_file(all_fasta,VH_folder,type_char,mkdir)
cp_all_file(hmm_out_folder,VH_folder,type_char,mkdir)
type_char="VL"
mkdir="F"
cp_all_file(all_fasta,VL_folder,type_char,mkdir)
cp_all_file(hmm_out_folder,VL_folder,type_char,mkdir)
#2.�ֱ��VH��VL���з��������ȷ���VH
scheme="Longest_CDR_H"
dict_match=match_fasta_txt(VH_folder)
print dict_match
get_cdr_from_hmmer(VH_folder,dict_match,scheme)
scheme="Longest_CDR_L"
dict_match=match_fasta_txt(VL_folder)
print dict_match
get_cdr_from_hmmer(VL_folder,dict_match,scheme)
#3.ȡ��CDR�����ļ����ֱ������VH_VL�Լ�VH�������ļ����У����ں�����ʹ��    
type_char="CDR_sequences"
mkdir="F"
cp_all_file(VH_folder,VH_VL_folder,type_char,mkdir)
cp_all_file(VL_folder,VH_VL_folder,type_char,mkdir)
cp_all_file(VH_folder,VH_CDR3,type_char,mkdir)
#################################################���Ĳ�########################################################
#######����merge_VH_VL��VH��VL�����е�CDR����ϲ�
#1.���Ƚ����������������п�����combined_folder��ȥ�����ں����ĺϲ�����Ľ���
pos_char=".fasta"#��ȫ����fasta�ļ��������µ��ļ�����ȥ
pos="end"
cp_start_end(VH_VL_folder,combined_folder,pos_char,pos)
pos_char="C"#���ｫC��ͷ���ļ�������δ����Ҫȥ�����ֵ��ļ�����ȥ
pos="start"
cp_start_end(VH_VL_folder,C_combined,pos_char,pos)

#2.Ȼ����кϲ���VH+VL���ȴ������е����ݣ��ڴ���C��ͷ������
#���Ƚ�C��ͷ���ļ����д��������Ĳ����ͳ����һ�������ڳ���ĺϲ�
species="rabbit"
divide_species_cdr(C_combined,species)
type_char=("without")
mkdir="F"
cp_all_file(C_combined,rm_folder,type_char,mkdir)
#���ļ���ID���д����Ա�֤�ںϲ���ʱ�򲻻��������
rm_file="T"
split_char="_"
pos_start=0
pos_end=-3
split_id(combined_folder,split_char,pos_start,pos_end,rm_file)
merge_VH_VL(combined_folder)
split_id(rm_folder,split_char,pos_start,pos_end,rm_file)
merge_VH_VL(rm_folder)
#�½��ļ��д���C��ͷ��VH_VL�ϲ��ļ��Լ����кϲ��ļ�
type_char="combined_VH_VL_CDR"
mkdir="F"
cp_all_file(rm_folder,C_combined_all,type_char,mkdir)
cp_all_file(combined_folder,combined_all,type_char,mkdir)

############################################���岽#########################################
###############################################################################################
#########################������������CDR3��������ȡ����#################################
#1.��������C��ͷ��CDR3���ļ��У����ҿ�����ȥC��ͷ��CDR3

#2.��C��ͷ����������ת�Ƶ�C_CDR3��
type_char="C"
mkdir="F"
cp_all_file(VH_CDR3,C_CDR3,type_char,mkdir)
#3.�������ֽ��з���
species="rabbit"
divide_species_cdr(C_CDR3,species)
#4.�½��ļ��д��治���ض����ֵ�����CDR3
type_char="without"
mkdir="F"
cp_all_file(C_CDR3,rm_folder_CDR3,type_char,mkdir)
#��ȡCDR3,���Ҵ�����CDR3_only�ļ����У������ܵ�CDR3��ֻ��C��ͷ�ģ�
get_cdr3(rm_folder_CDR3)
get_cdr3(VH_CDR3)
type_char="CDR3_only"
mkdir="F"
cp_all_file(rm_folder_CDR3,species_CDR3_only_folder,type_char,mkdir)
cp_all_file(VH_CDR3,CDR3_only_folder,type_char,mkdir)
#�༭ID
rm_file="T"
split_char="_"
pos_start=0
pos_end=-3
split_id(species_CDR3_only_folder,split_char,pos_start,pos_end,rm_file)
split_id(CDR3_only_folder,split_char,pos_start,pos_end,rm_file)
#�����༭ID������п������µ��ļ�����
type_char="edited_id"
mkdir="F"
cp_all_file(species_CDR3_only_folder,species_CDR3_formal_id_folder,type_char,mkdir)
cp_all_file(CDR3_only_folder,CDR3_formal_id_folder,type_char,mkdir)
#�����е���Ҫ����Ƶ�ʷ������ļ��н���Ƶ�ʷ���
#��һ�����ǶԳ������еķ���
all_path=[species_CDR3_formal_id_folder,CDR3_formal_id_folder,C_combined_all,combined_all]
for path_n in all_path:
    if path_n.find("CDR3")!=-1:
        len_out=1
        freq_analysis_fasta_ID(path_n,len_out)
    elif path_n.find("all_combined")!=-1:
        len_out=6
        freq_analysis_fasta_ID(path_n,len_out)
#����Ӧ���ļ���������Ӧ���ļ�����
cp_file_to_result(file_repo)
