#coding=utf-8
from Bio.SearchIO.HmmerIO import Hmmer3TextParser
import os
from Bio import SeqIO
from os import path,access,R_OK
##########
#PART##
###########
#这一部分主要是用于根据给出的条件运行HMMER












##########
#PART2##
###########
#这一部分主要是解析HMMER的输出文件

#function1:解析hmmer输出文件
def hmmer_parse(path,hmmer_output):
    os.chdir(path)
    file_iter=open(hmmer_output,"r")#把输出文件读成一个句柄
    file_n=Hmmer3TextParser(file_iter)#利用HmmerTextParser对句柄进行解析
    return file_n

#function2：解析fasta文件
#缓存序列文件,把其储存在字典中，字典的键为序列的ID，值为碱基序列
def read_fasta(fasta_file):
    sequences={}#首先定义一个空的字典，用于后续储存序列
    all_seq=SeqIO.parse(fasta_file,"fasta")#利用SeqIO对fasta文件进行读取
    #遍历fasta文件，并且将序列储存在字典中
    for seq in all_seq:
        seq_id=str(seq.id)
        seq_seq=str(seq.seq)
        sequences[seq_id]=seq_seq
    #返回保存有序列的字典
    return sequences
#function3:对解析后的hmmer柄文件进行分析
def hmmer_handle_analysis(hmmer_parse):
    #首先检查比对文件中有几个domain,再分别对domain进行分析
    #print annotations
    anno_seq=[]#为每个序列的注释定义一个列表
    anno_domain=[]#为每一个序列所有的结构域定义一个列表，这个列表将来包含在anno_seq中
    
    for annotations in hmmer_parse:
        print annotations
        anno_domain_part=[]#这个列表将会记录单个结构域的注释信息，这个列表将会包含在anno_domain中
        for query in annotations:#这里的query也只是代表一个序列
            #print query
            for hsp in query.hsps[0]:#对每一个hsp进行遍历,记录其中的详细信息，hmm开始，结束，query开始结束
                query_id=hsp.query_id
                hmm_start=hsp.hit_start
                hmm_end=hsp.hit_end
                query_start=hsp.query_start
                query_end=hsp.query_end        
                hmm_stat_origin=hsp.aln_annotation["RF"]#取出来hmm状态的标志序列
                #print hmm_stat_origin
                query_stat_origin=hsp.aln_annotation["PP"]#取出query序列状态的标志序列
                len_insert=0
                len_delete=0
                insert_pos_list=[]#记录所有的插入的连续的位置
                insert_pos=[]#记录插入起始的位置
                insert_len=[]#记录插入序列的长度
                all_insert=[]#将所有的插入序列都放入到一个列表中，用元组的形式储存信息
                del_pos_list=[]#记录所有的缺失的连续的位置
                del_pos=[]#记录缺失起始的位置
                del_len=[]#记录缺失序列的长度
                all_del=[]#将所有的缺失（起始和长度）都放入到一个元组中，并且都放置在这个列表中
                #首先利用一个for循环记录下来每个插入的序列的位置和长度
                for pos,stat in enumerate(hmm_stat_origin):
                    if stat==".":
                        insert_pos_list.append(pos+hmm_start)
                    elif stat=="x" and hmm_stat_origin[pos-1]==".":
                        insert_pos.append(insert_pos_list[0])
                        insert_len.append(insert_pos_list[-1]-insert_pos_list[0]+1)
                        insert_pos_list=[]
                    elif stat=="x" and hmm_stat_origin[pos-1]!=".":
                        insert_pos_list=[]
                #在利用一个for循环记录下来每个缺失位置的起始位置和长度
                for pos,stat in enumerate(query_stat_origin):
                    if stat==".":
                        del_pos_list.append(pos+hmm_start)
                    elif stat=="x" and query_stat_origin[pos-1]==".":
                        del_pos.append(del_pos_list[0])
                        del_len.append(del_pos_list[-1]-del_pos_list[0]+1)
                        del_pos_list=[]
                    elif stat=="x" and query_stat_origin[pos-1]!=".":
                        del_pos_list=[]
                if insert_pos and insert_len:
                    for i in xrange(len(insert_pos)):
                        all_insert.append((insert_pos[i],insert_len[i]))#利用元组的方式返回插入序列的开始位置和长度，记录在anno_domain中
                if del_pos and del_len:
                    for i in xrange(len(insert_pos)):
                        all_del.append((del_pos[i],del_len[i]))#利用元组的方式返回插入序列的开始位置和长度，记录在anno_domain中
                anno_domain_part.append((all_insert,all_del))#同时记录下每个结构域的开始和结束的位置以及序列匹配的开始和结束的位置
                #print anno_domain_part
                #anno_domain_part.append(all_del)
                #print anno_domain_part
        anno_domain.append((query_id,hmm_start,hmm_end,query_start,query_end,anno_domain_part))#最后一项为插入缺失位置和长度的记录，假如结构域为两个那么最后一项为多个列表的元组
        #print anno_domain
        #print "this is %s"%anno_domain
    anno_seq=anno_domain
    return anno_seq

#function4:结合hmmer柄文件分析的结果,得出不同结构域的起始和结束的位置

def get_pos_domain(anno_seq,seq_dict,scheme="IMGT"):#目前只考虑一个domain的情况
    schemes_limit={"IMGT":[0,25,33,55,65,104,117],
              "AHO":[],
              "Kabat":[],
              "Chothia":[]}
    scheme_limit=schemes_limit[scheme.upper()]
    region_anti=["FR1","CDR1","FR2","CDR2","FR3","CDR3"]
    region_seq={}#用来记录每个序列的各个结构域的序列
    for seq_anno in anno_seq:
        hmm_start=seq_anno[1]#hmm_start,hmm_end,query_start,query_end在目前仅考虑只有一个结构域的情况
        hmm_end=seq_anno[2]
        query_start=seq_anno[3]
        query_end=seq_anno[4]
        if len(seq_anno[-1])>1:#取出含有插入和缺失位置信息的列表，列表内容为'M03186:269:000000000-BK7T5:1:1101:17338:2018:CCTCCTGA+TCTTTCCC', 0, 112, 0, 121, [([(98, 9)], [])]的最后一项
            domain_details=seq_anno[-1][0]
        elif len(seq_anno[-1])==1:
            domain_details=seq_anno[-1]
        domain_num=1#目前仅遍历第一个结构域
        insert_num=len(domain_details[0][0])
        del_num=len(domain_details[0][1])#domain[0代表取出第一个domain]
        #domain_num=len(domain_details)
        len_to_add=min(hmm_start,query_start)#确定最小的开始位置，以免后续获取序列的时候
        if hmm_start<10:
            seq_start=query_start-hmm_start#找到序列开始位置的索引，这也是FR1开始的位置
        elif hmm_start>10:
            pass#暂时还不知道如何处理
        
        query_id=seq_anno[0]
        #print seq_anno
        insertion=[]
        deletion=[]
        #取出各个部分的序列，放在query_seq中
        print query_id
        query_seq=seq_dict[query_id]
        print query_seq
        #利用domain_number对所有的结构域进行分析
        print domain_details
        for num in xrange(domain_num):#对所有的domain进行遍历提取出所有的插入缺失的位置，目前仅对一个结构域进行遍历
            #判断每个插入是不是为零，假如不为零，那么就遍历取出来，假如为零的pass
            if len(domain_details[num][0])>0:
                for insert_x in domain_details[num][0]:
                    insertion.append(insert_x)#获取当前所有的插入位置和长度
            if len(domain_details[num][1])>0:
                for delete_x in domain_details[num][1]:
                    deletion.append(delete_x)#获取当前所有位置的缺失位置和长度
        domain_seq={}
        print insertion,deletion
        for i in xrange(len(scheme_limit)-1):#遍历每一个区间,减一是因为最后一个是 CDR3的末端，没有成对的匹配
            start=scheme_limit[i]#开始的位置是第一个
            end=scheme_limit[i+1]#结束的位置是开始的下一个位置
            len_insertion=0#len_insertion是用来记录在结束的位置之前插入碱基的长度
            len_deletion=0#len_deletion是用来记录在结束的位置之前缺失碱基的长度
            if len(insertion)>0 and len(deletion)>0:#假如该序列同时有插入跟缺失
                for insert_n in insertion:#利用for循环累计在结束位置之前总的插入长度和缺失长度
                    if insert_n[0]<end:
                        len_insertion+=insert_n[1]
                    elif insert_n>=end:
                        pass
                for del_n in deletion:
                    if del_n[0]<end:
                        len_deletion+=del_n[1]
                    elif del_n>=end:
                        pass
            
            elif len(insertion)>0 and len(deletion)==0:
                for insert_n in insertion:#利用for循环累计在结束位置之前总的插入长度和缺失长度
                    print"ooooooooooooooook"
                    print insert_n
                    if insert_n[0]<end:
                        len_insertion+=insert_n[1]
                        print "ooooooooooooooo%d"%len_insertion
                    elif insert_n>=end:
                        pass


            elif len(insertion)==0 and len(deletion)>0:
                for del_n in deletion:
                    if del_n[0]<end:
                        len_deletion+=del_n[1]
                    elif del_n>=end:
                        pass
            

            elif len(insertion)==0 and len(deletion)==0:
                pass
            print len_insertion,len_deletion
            if i ==0:
                start_pos=query_start-len_to_add#获取序列的起始位置
            elif i>0:
                start_pos=end_pos
            
            end_pos=query_start+(end-hmm_start)+len_insertion-len_deletion
            print start_pos
            print end_pos
            print type(start_pos)
            print type(end_pos)
            sequence=query_seq[start_pos:end_pos]
            domain_seq[region_anti[i]]=sequence
            print domain_seq
            len_insertion=0
            len_deletion=0
        region_seq[query_id]=(domain_seq)
    return region_seq
        
            
##########
#PART3##
###########
#这一部分主要是对获取的序列进行编号               
            
            
            
        
            







#首先第一步解析hmm输出文件，并且缓存序列文件
#首先转入到目标文件夹，并且读取相应的hmmer输出的文件
#下面是自定义文件名和路径名，file_path,hmmer_output,fasta_file分别是hmmer输出文件所在的文件夹，hmmer输出的文件名，含有氨基酸序列的fasta文件名
if __name__=="__main__":
    file_path="E:\HMM_transfer\HMM_transfer_3ed"
    hmmer_output="VH_output_little.txt"
    fasta_file="Genesp-VH_AA.fasta"
    #利用判断语句检测目标文件夹目录中有没有该文件以及fasta_file是不是fasta文件
    hmmer_file=path.join(file_path,hmmer_output)
    fasta_name,postfix=os.path.splitext(fasta_file)
    assert path.exists(file_path) and path.exists(hmmer_file),"there are something wrong with the hmmer_output directory and file"
    assert postfix==".fasta","the format of file should be fasta"

    #读取hmmer的输出文件，并且返回hmmer输出文件的解析文件
    hmmer_handle=hmmer_parse(file_path,hmmer_output)
    #缓存序列文件
    seq_dict=read_fasta(fasta_file)
    #利用hmmer输出文件的解析和fasta文件进行分析
    annotation_all_seq=hmmer_handle_analysis(hmmer_handle)
    print annotation_all_seq
    domain_sequences=get_pos_domain(annotation_all_seq,seq_dict)


   

