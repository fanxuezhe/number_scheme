#coding=utf-8
from Bio.SearchIO.HmmerIO import Hmmer3TextParser
import os
from Bio import SeqIO
from os import path,access,R_OK

##########################
#####FGS_ANALYSIS#########
##########################

#定义一个function,该function的功能是通过给定路径和引物（或者是需要的序列的前面的序列）以及VH，VL之间的linker序列，提取一代测序结果中的核苷酸序列，氨基酸序列，并且根据将其分成VH,VL两个序列，并且分别保存在不同的文件夹中
#输入的参数：路径，引物（正向和反向），linker序列
#输出的结果：保存在DNA_seq，AA_seq,VL_seq,VH_seq中的分别是带有引物的核苷酸序列，氨基酸序列的全长，VL和VH序列
def FGS_analysis(path,primer_f,primer_r):
    import os
    from Bio.Seq import Seq
    os.chdir(file_path)
    path_name=file_path.split(os.sep)[-1]
    DNA_seq=open("%s_sequences_DNA.fasta"%path_name,"w")
    AA_seq=open("%s_sequences_AA.fasta"%path_name,"w")
    VL_seq=open("%s_sequences_VL.fasta"%path_name,"w")
    VH_seq=open("%s_sequences_VH.fasta"%path_name,"w")
    all_files=os.listdir(file_path)
    i=1
    print all_files
    for file_n in all_files:
        file_name,file_type=os.path.splitext(file_n)
       
        #print file_type
        #print file_name
        vl_seq=""#将循环值归零，记得假如要判断循环中的值是否为真，那么该值一定要在循环中归零
        vh_seq=""
        if file_type==".seq":
            print file_n
            file_content=open(file_n,"r").readlines()
            #由于seq文件中的序列是多行分割的因此，需要把序列进行拼接
            all_seq=[]
            for seq_n in file_content:
                seq_short=seq_n.strip("\n ")
                all_seq.append(seq_short)
            sequence_all="".join(all_seq)
            print sequence_all
            aa_seq=""
            #print sequence_all
            print file_name
            for f_n in primer_f:
                if sequence_all.find(f_n)!=-1:
                    sequence_all_1=sequence_all.split(f_n)[1]
                    #print sequence_all_1
                    for r_n in primer_r:
                        print r_n
                        if sequence_all_1.find(r_n)!=-1:
                            aa_seq=str(Seq(sequence_all_1.split(r_n)[0]).translate())
                        
                            break
            assert aa_seq,"there are primers that need to be added"
            if aa_seq.find("*")==-1:
            #print sequence_all
                print aa_seq
                print file_name
                #将DNA序列写到fasta文件中
                DNA_seq.write(">"+str(file_name)+str(i)+"\n")
                DNA_seq.write(sequence_all+"\n")
                #将氨基酸序列写到fasta文件中
                AA_seq.write(">"+str(file_name)+str(i)+"\n")
                AA_seq.write(aa_seq+"\n")
                #将氨基酸序列分为两个部分，VL和VH
                for sp_n in split_all:
                    if aa_seq.find(sp_n)!=-1:
                        vl_seq,vh_seq=aa_seq.split(sp_n)
                        name_id="".join(file_name.split(" "))#由于当id中存在空格的时候，hmmscan处理数据时会将id自动断开，取第一部分作为id,因此去除掉
                        print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX%d"%i
                        print vl_seq
                        print vh_seq
                        VL_seq.write(">"+name_id+str(i)+"\n")
                        VL_seq.write(vl_seq+"\n")
                        VH_seq.write(">"+name_id+str(i)+"\n")
                        VH_seq.write(vh_seq+"\n")
                        break
                print vl_seq
                print vh_seq
                assert vl_seq and vh_seq,"there is linker need to add"
                
                #这里利用i来记录这是第几条序列，同时也用来表示fasta文件中序列ID的唯一性
                i+=1
                print i
            else:
                print "there is sequence contain the *"
            
        
    DNA_seq.close()           
    AA_seq.close()
    VH_seq.close()
    VL_seq.close()
#定义一个新的功能，这个程序能够记录所有含有终止密码子的氨基酸序列，这些氨基酸序列可能是因为测序的原因导致最后的终止密码子的出现
#输入的参数：路径，引物（正向和反向）
def record_abnormal(path,primer_f,primer_r):
    import os
    from Bio.Seq import Seq
    os.chdir(path)
    files_all=os.listdir(path)
    file_abnormal=open("all_abnormal_sequences.fasta","w")
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        #判断要处理的文件是不是以seq结尾的文件
        if file_type==".seq":
            #接下来合并所有的核苷酸序列，然后翻译成氨基酸
            file_content=open(file_n,"r").readlines()
            sequence="".join([seq_n.strip("\n") for seq_n in file_content])
            print sequence
            aa_seq=""
            for f_n in primer_f:#遍历每一个引物，假如引物找不到那么就记录下来
                if sequence.find(f_n)!=-1:
                    sequence_all_1=sequence.split(f_n)[1]
                    #print sequence_all_1
                    for r_n in primer_r:
                        print r_n
                        if sequence_all_1.find(r_n)!=-1:
                            aa_seq=str(Seq(sequence_all_1.split(r_n)[0]).translate())
                        
                            break
            assert aa_seq,"there are primers that need to be added"#检测是否正常的检测出序列，假如检测不出来那么就说明引物有问题，需要添加新的引物
            if aa_seq.find("*")!=-1:
                file_abnormal.write(file_n+"\n")
                file_abnormal.write(aa_seq+"\n")
    file_abnormal.close()



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
#接下来定义一个新的功能，该功能主要是利用给定的路径,将其中的fasta文件（即氨基酸序列文件）和hmmer的输出文件进行配对
def  match_fasta_txt(file_path):   
    os.chdir(file_path)
    files_all=os.listdir(file_path)
    dict_txt_fasta={}
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        fasta_name=file_name+"."+"fasta"
        print fasta_name
        print file_name
        print file_type
        if file_type==".txt" and files_all.count(fasta_name)==1:
            dict_txt_fasta[file_n]=fasta_name#将配对成功的文件储存在字典中，方便以后调用
        print dict_txt_fasta#输出能够成功配对的fasta和Txt文件

    return dict_txt_fasta
#接下来定义一个功能，该功能利用match_fasta_txt功能返回的字典（含有序列和hmmer输出文件名称的配对）进行分析，获得含有cdr序列文件,"""需要输入的参数是：文件所在的文件夹和返回的字典数据以及要选择的scheme种类"""
def get_cdr_from_hmmer(file_path,dict_txt_fasta,scheme="IMGT"):
    os.chdir(file_path)
    for txt_file,fasta_file in dict_txt_fasta.items():#对配对的fasta和txt文件进行遍历
        file_name=os.path.splitext(txt_file)[0]#获得txt文件的文件名
        file_cdr_seq=open("%s_CDR_sequences.fasta"%file_name,"w")#利用txt文件的文件名明明新的文件，它将会储存cdr序列的fasta文件
        hmmer_output=txt_file#定义txt文件名，用于后续分析
        fasta_file=fasta_file#定义fasta文件名用于后续分析
        #利用判断语句检测目标文件夹目录中有没有该文件以及fasta_file是不是fasta文件
        hmmer_file=os.path.join(file_path,hmmer_output)
        fasta_name,postfix=os.path.splitext(fasta_file)
        assert os.path.exists(file_path) and os.path.exists(hmmer_file),"there are something wrong with the hmmer_output directory and file"#利用assert语句进行判断hmmer输出文件是否存在
        assert postfix==".fasta","the format of file should be fasta"#检查fasta文件的后缀是不是fasta

        #读取hmmer的输出文件，并且返回hmmer输出文件的解析文件
        hmmer_handle=hmmer_parse(file_path,hmmer_output)#都没问题以后
        #缓存序列文件
        seq_dict=read_fasta(fasta_file)
        #利用hmmer输出文件的解析和fasta文件进行分析
        annotation_all_seq=hmmer_handle_analysis(hmmer_handle)
        print annotation_all_seq
        domain_sequences=get_pos_domain(annotation_all_seq,seq_dict,scheme)
        print domain_sequences
        for query_id,domain_seq in domain_sequences.items():#获取对应的CDR序列
            CDR1_seq=domain_seq["CDR1"]
            CDR2_seq=domain_seq["CDR2"]
            CDR3_seq=domain_seq["CDR3"]
            file_cdr_seq.write(">"+query_id+"\n")
            file_cdr_seq.write(CDR1_seq+"\t"+CDR2_seq+"\t"+CDR3_seq+"\n")#将CDR序列写入文件中
        file_cdr_seq.close()




            
#function3:对解析后的hmmer柄文件进行分析
def hmmer_handle_analysis(hmmer_parse):
    #首先检查比对文件中有几个domain,再分别对domain进行分析
    #print annotations
    anno_seq=[]#为每个序列的注释定义一个列表
    anno_domain=[]#为每一个序列所有的结构域定义一个列表，这个列表将来包含在anno_seq中
    
    for annotations in hmmer_parse:
        #print dir(annotations.hsps)
        anno_domain_part=[]#这个列表将会记录单个结构域的注释信息，这个列表将会包含在anno_domain中
        #for query in annotations:#这里的query也只是代表一个序列
        #print query
        i=0#控制仅取第一个最好的匹配
        for hsp in annotations.hsps:#对每一个hsp进行遍历,记录其中的详细信息，hmm开始，结束，query开始结束
            if i<1:
                #print i
                query_id=hsp.query_id
                #print "query_id%s"%(query_id)
                hmm_start=hsp.hit_start
                hmm_end=hsp.hit_end
                query_start=hsp.query_start
                query_end=hsp.query_end        
                hmm_stat_origin=hsp.aln_annotation["RF"]#取出来hmm状态的标志序列
                print hmm_stat_origin
                query_stat_origin=hsp.aln_annotation["PP"]#取出query序列状态的标志序列
                print query_stat_origin
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
                        #print "kkkkkkkkkkkkkkkkkkkkkk"
                        del_pos_list.append(pos+hmm_start)
                    elif stat!="." and query_stat_origin[pos-1]==".":
                        del_pos.append(del_pos_list[0])
                        del_len.append(del_pos_list[-1]-del_pos_list[0]+1)
                        del_pos_list=[]
                    elif stat!="." and query_stat_origin[pos-1]!=".":
                        del_pos_list=[]
                if insert_pos and insert_len:
                    for i in xrange(len(insert_pos)):
                        all_insert.append((insert_pos[i],insert_len[i]))#利用元组的方式返回插入序列的开始位置和长度，记录在anno_domain中
                if del_pos and del_len:
                    for i in xrange(len(del_pos)):
                        all_del.append((del_pos[i],del_len[i]))#利用元组的方式返回插入序列的开始位置和长度，记录在anno_domain中
                anno_domain_part.append((all_insert,all_del))#同时记录下每个结构域的开始和结束的位置以及序列匹配的开始和结束的位置
                i+=1
                print anno_domain_part
                #anno_domain_part.append(all_del)
                #print anno_domain_part
        anno_domain.append((query_id,hmm_start,hmm_end,query_start,query_end,anno_domain_part))#最后一项为插入缺失位置和长度的记录，假如结构域为两个那么最后一项为多个列表的元组
        #print anno_domain
        #print "this is %s"%anno_domain
    anno_seq=anno_domain
    return anno_seq

#function4:结合hmmer柄文件分析的结果,得出不同结构域的起始和结束的位置

def get_pos_domain(anno_seq,seq_dict,scheme="IMGT"):#目前只考虑一个domain的情况
    schemes_limit={"IMGT":[0,26,38,55,65,104,117],
              "AHO":[],
              "Kabat":[],
              "Chothia":[],
            "LONGEST_CDR_H":[0,26,40,54,74,104,117],
            "LONGEST_CDR_L":[0,23,40,55,69,104,117]     }
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
def merge_VH_VL(file_path):#这部分功能主要是将VH和VL的CDR区进行合并随后在进行一个频率分析
    os.chdir(file_path)
    files_all=os.listdir(file_path)
    dict_cdr_vh_vl={}
    for file_n in files_all:
         if file_n.find("CDR_sequences")!=-1 and file_n.find("VL")!=-1:#首先找到储存CDR序列的文件
             print file_n
             file_name=file_n.split("VL")[0]
             file_n_vh=file_n.replace("VL","VH")#这里定义vh文件的文件名称
             
             if files_all.count(file_n_vh)==1:
                 dict_cdr_vh_vl[file_n]=file_n_vh#先将成对的文件储存在字典中
                 print dict_cdr_vh_vl
    for vl_cdr,vh_cdr in dict_cdr_vh_vl.items():#对成对的文件进行分析
        vl_all=open(vl_cdr,"r").readlines()#读取内容
        vh_all=open(vh_cdr,"r").readlines()
        assert len(vl_all)==len(vh_all),"attention the number of CDR of VL is different from the numner of CDR of VH"
        dict_vl={}
        dict_vh={}
        for i in range(0,len(vl_all),2):#这个循环的目的是把VH和VL的id和序列都储存在字典中，这样就可以方便后续用ID进行匹配序列
            print vl_all[i]
            print vl_all[i+1]
            if vl_all[i][0]==">" and vl_all[i+1][0]!=">":
                dict_vl[vl_all[i]]=vl_all[i+1]
            if vh_all[i][0]==">" and vh_all[i+1][0]!=">":
                dict_vh[vh_all[i]]=vh_all[i+1]
            else:
                print "the fasta file is not standard the sequence not follow the id"
        #接下来要合并两个文件
        file_name=os.path.splitext(vl_cdr)[0]
        combine_cdr=open("%s_combined_VH_VL_CDR.fasta"%file_name,"w")
        for vl_id,vl_seq in dict_vl.items():
            combined_seq=vl_seq.strip("\n")+"\t"+dict_vh[vl_id]
            combine_cdr.write(vl_id)
            combine_cdr.write(combined_seq)
        combine_cdr.close()

#接下来写一个function，用来统计fasta文件中序列的频率，该程序将会统计所有的fasta文件的频率
def freq_analysis_fasta(path):
    import os,shutil
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_type=os.path.splitext(file_n)[1]
        if file_type==".fasta":
            all_seq=open(file_n,"r").readlines()
            unique_seq=[]
            dict_unique_seq={}#利用一个字典对唯一序列进行记录这样方便后续对ID的书写
            freq_seq={}
            header=file_n.split(".")[0]
            file_freq=open("%s_freq.fasta"%header,"w")
            file_freq.write("ID"+"\t"+"freq_of_seq"+"\t"+"sequences"+"\n")
            #记录所有的唯一的DNA序列
            for i in range(1,len(all_seq),2):
                print all_seq[i]
                if unique_seq.count(all_seq[i])==0:
                    unique_seq.append(all_seq[i])
                    dict_unique_seq[all_seq[i-1]]=all_seq[i]#将ID和序列记录下来,id含有回车符号
            for id_n,seq_n in dict_unique_seq.items():
                seq_freq=all_seq.count(seq_n)
                seq=seq_n.strip("\n")#这里的序列含有回车符号
                freq_seq[id_n]=seq_freq
            sorted_seq=sorted(freq_seq.items(),key=lambda item:item[1],reverse=True)#这个将按照字典的值的大小进行排序,返回的sorted_seq是一个列表类型
            print sorted_seq
            for id_n,freq_n in sorted_seq:#遍历过程中取出id,频率
                seq_n=dict_unique_seq[id_n]
                file_freq.write(id_n.strip("\n")+"\t"+str(freq_n)+"\t"+seq_n.strip("\n")+"\n")
            file_freq.close()

#重新定义一个function这个用来分析频率以后能够将该频率下的所有ID分别记下来
def freq_analysis_fasta_ID(path):
    import os,shutil
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_type=os.path.splitext(file_n)[1]
        if file_type==".fasta":
            all_seq=open(file_n,"r").readlines()
            unique_seq=[]
            dict_unique_seq={}#利用一个字典对唯一序列进行记录这样方便后续对ID的书写
            freq_seq={}
            header=file_n.split(".")[0]
            file_freq=open("%s_freq.fasta"%header,"w")
            file_freq.write("freq_of_seq"+"\t"+"sequences"+"\t"+"IDs"+"\n")
            #记录所有的唯一的DNA序列
            for i in range(1,len(all_seq),2):
                print all_seq[i]
                if unique_seq.count(all_seq[i])==0:
                    unique_seq.append(all_seq[i])
                    dict_unique_seq[all_seq[i-1]]=all_seq[i]#将ID和序列记录下来,id含有回车符号
            for id_n,seq_n in dict_unique_seq.items():
                ids=[]
                for m in xrange(len(all_seq)):
                    if all_seq[m]==seq_n:
                        ids.append(all_seq[m-1])
                seq_freq=all_seq.count(seq_n)
                seq=seq_n.strip("\n")#这里的序列含有回车符号
                freq_seq[id_n]=(seq_freq,ids)
            sorted_seq=sorted(freq_seq.items(),key=lambda item:item[1][0],reverse=True)#这个将按照字典的值的大小进行排序,返回的sorted_seq是一个列表类型
            print sorted_seq
            for id_n,freq_n in sorted_seq:#遍历过程中取出id,频率
                seq_n=dict_unique_seq[id_n]
                freq=freq_n[0]
                ids=[x.strip("\n")  for x in freq_n[1]]
                file_freq.write(str(freq)+"\t"+str(seq_n).strip("\n")+"\t"+"  ".join(ids)+"\n")
            file_freq.close()
#该程序的功能是针对特定序列有部分杂序列的清除,header_str代表的是一个分隔符，分隔符前面的序列以及分隔符将会被删除
def remove_header_vh(file_path,header_str):
    os.chdir(file_path)
    files_all=os.listdir(file_path)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta" and file_name.endswith("VH"):
            file_content=open(file_n,"r").readlines()
            file_edited=open("%s_edited.fasta"%file_name,"w")#去掉接头的序列后新建文件储存序列并且在文件名后方加上edited
            for seq_n in file_content:
                if not seq_n.startswith(">") and seq_n.find(header_str)!=-1:
                    seq_edited=seq_n.split(header_str)[1]#这里用1表示只考虑了去掉N端的多余序列
                    file_edited.write(seq_edited)#由于没有处理后边的回车符，这里就不在加上去了。
                elif seq_n.startswith(">"):
                    file_edited.write(seq_n)
            file_edited.close()
                    

#接下来定义一个新的功能，该功能用来将fasta文件的id缩短，目前的程序以“-”为分割，并且取前四个元素作为未来的id,split_char代表的是需要使用到的分隔符，pos_start,pos_end代表的是分隔符开始和结束的位置，rm_file代表的是是否要移除原来的文件,用T来表示即可,默认是不删除的
#注意到，假如要取最后一位的时候需要用[-1:None]
def split_id(path,split_char,pos_start,pos_end,rm_file="F"):
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta":
            file_content=open(file_n,"r").readlines()
            new_file=open("%s_edited_id.fasta"%file_name,"w")
            for line_n in file_content:
                if line_n.startswith(">"):
                    new_line=split_char.join(line_n.split(split_char)[pos_start:pos_end])
                    if new_line.startswith(">"):#这里判断剪切之后是否还是以">"结尾，假如不是">"结尾，那么就需要重新加上去
                        pass
                    elif not new_line.startswith(">"):
                        new_line=">"+str(new_line)#要确保两者都是字符，以免相加的时候出现问题
                    new_file.write(new_line+"\n")
                else:
                    new_file.write(line_n)#由于该行本身含有换行符
            new_file.close()
            if rm_file=="T":
                os.remove(file_n)
        
#这里定义一个模块用来构建氨基酸序列的发育树
def phylo_draw(path):
    from align_own import clustalw_align_alone as clu
    from Bio import Phylo as pl
    import pylab
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta":
            clu(path,file_n)
            dnd_file=file_name+".dnd"
            print dnd_file
            tree=pl.read(dnd_file,"newick")
            #tree.root=True
            #pl.draw_graphviz(tree,prog="dot")
            pl.draw(tree,branch_labels=lambda c:c.branch_length )
            pylab.savefig("%s_fig.png"%file_name)
#定义一个新功能，该功能的目的是删除掉fasta文件中重复的序列，获得的唯一序列用来做进化树分析
def rm_repeat_seq(path,rm_file="F"):
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta":
            file_content=open(file_n,"r").readlines()
            file_new=open("%s_rm_repeat.fasta"%file_name,"w")
            unique_seq=[]#记录唯一的序列
            for i in xrange(len(file_content)):
                line_n=file_content[i]
                if not line_n.startswith(">") and unique_seq.count(line_n)==0:
                        unique_seq.append(line_n)
                        file_new.write(file_content[i-1])
                        file_new.write(line_n)
            file_new.close()
            if rm_file=="T":#检测是否需要删除源文件
                os.remove(file_n)
#接下来定义一个新的function，该功能是利用已有的含有CDR1/2/3的fasta文件转化成为含有CDR3的序列,序列的ID不变
def get_cdr3(path,rm_file="F"):
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta":
            file_content=open(file_n,"r").readlines()
            file_new=open("%s_CDR3_only.fasta"%file_name,"w")
            for i in xrange(len(file_content)):
                line_n=file_content[i]
                if line_n.startswith(">"):
                    file_new.write(line_n)
                elif not line_n.startswith(">"):
                    file_new.write(line_n.split("\t")[-1])
            file_new.close()
            if rm_file=="T":#检测是否需要删除源文件
                os.remove(file_n)
#定义一个function，该function的功能是合并多个fasta文件，
def combined_fasta(path):
    os.chdir(path)
    files_all=os.listdir(path)
    file_new=open("conbined_all_seq.fasta","w")
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta":
            file_content=open(file_n,"r").readlines()
            for seq_n in file_content:
                file_new.write(seq_n)
        
    file_new.close()


#定义一个function,该function的功能是统计频率后将频率作为id，将序列作为序列
def freq_analysis_num_id(path):
    import os,shutil
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_type=os.path.splitext(file_n)[1]
        if file_type==".fasta":
            all_seq=open(file_n,"r").readlines()
            unique_seq=[]
            dict_unique_seq={}#利用一个字典对唯一序列进行记录这样方便后续对ID的书写
            freq_seq={}
            header=file_n.split(".")[0]
            file_freq=open("%s_freq.fasta"%header,"w")
            #记录所有的唯一的DNA序列
            for i in range(1,len(all_seq),2):
                print all_seq[i]
                if unique_seq.count(all_seq[i])==0:
                    unique_seq.append(all_seq[i])
                    dict_unique_seq[all_seq[i-1]]=all_seq[i]#将ID和序列记录下来,id含有回车符号
            #统计频率
            for id_n,seq_n in dict_unique_seq.items():
                seq_freq=all_seq.count(seq_n)
                seq=seq_n.strip("\n")#这里的序列含有回车符号
                freq_seq[id_n]=seq_freq
            sorted_seq=sorted(freq_seq.items(),key=lambda item:item[1],reverse=True)#这个将按照字典的值的大小进行排序,返回的sorted_seq是一个列表类型
            print sorted_seq
            i=1
            for id_n,freq_n in sorted_seq:#遍历过程中取出id,频率
                seq_n=dict_unique_seq[id_n]
                file_freq.write(">"+"num_"+str(i)+"\t"+"freq="+str(freq_n)+"\n")
                file_freq.write(seq_n.strip("\n")+"\n")
                i+=1
            file_freq.close()
#定义一个function，该function的功能是利用clustalw对文件夹中所有的氨基酸序列进行比对
def clustal_all_folder(path):
    from align_own import clustalw_align_alone
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_type=os.path.splitext(file_n)[1]
        if file_type==".fasta":
            clustalw_align_alone(path,file_n)
    
    
            


        
                 
        
            
##########
#PART3##
###########
#这一部分主要是对获取的序列进行编号               
        
            
            
        
            







#首先第一步解析hmm输出文件，并且缓存序列文件
#首先转入到目标文件夹，并且读取相应的hmmer输出的文件
#下面是自定义文件名和路径名，file_path,hmmer_output,fasta_file分别是hmmer输出文件所在的文件夹，hmmer输出的文件名，含有氨基酸序列的fasta文件名
if __name__=="__main__":
    import shutil,os
    #第三次运行
    #第一步利用FGS模块将核苷酸序列合并。并且拆分为VH VL 两个文件，同时记录下来含有终止密码子的序列
    #提前定义好特异性的参数
    #file_path=r"E:\Program\project\yanli\FGS_2018.11.8\051602-6 analysis"
    #F=["CAGTTTTGGCACAGGCGGCC","CAGAATTGGCCCAGGCGGCC"]
    #R=["GGGCCGTCGGTGGTC","GGGCCGTCGTGGTCG","GGCCGTCGTGGTCGA","GGCCGTCGGTGGTCG","GGGCCGTTGGTGGTCGACCT","GGCCCGTCGGTGGTCGACCT"]
    #split_all=["RAGGGGSGGGGSGGLGGAAA","RAGGGGSGGGGGGAAA","RAGGEGSGGGGGGAAA","APEEEGSGGGGSGGLGGAAA","RAGGGGSGGGGSGGLCGAAA","RAGGGGSEGGGSGGLGGAAA","RAGEGGSGGGGSGGLGGAAA","RAEGRGSGGGGSGGLGGAAA","RAGGEGSGGGGSGGLGGAAA"]#这里需要在循环过程中不断地添加新的linker序列


    #正式进行分析
    #FGS_analysis(file_path,F,R)
    #record_abnormal(file_path,F,R)
    #file_path=r"E:\Program\project\yanli\FGS_2018.11.8\20181029seq"
    #FGS_analysis(file_path,F,R)
    #record_abnormal(file_path,F,R)

    #第二步利用得到的VH和VL序列在linux系统中利用hmmscan进行分析


    #第三步利用match_fasta_txt匹配路径中的fasta和txt的匹配，并且利用匹配到的结果进行分析，得到CDR区域的序列
    ######################VH
    #path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\VH"
    #scheme="Longest_CDR_H"
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)
    #######################VL
    #path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\VL"
    #scheme="Longest_CDR_L"
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)
    #为了下一步的进行需要新建文件夹并且将所有的VH和VL的CDR合并
    #os.chdir(path)
    #path_to_mk=os.path.split(os.sep)
    #os.mkdir("")
    #第四步利用merge_VH_VL将VH和VL的序列的CDR区域合并
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\conbined_VH_VL"
    #merge_VH_VL(all_path)

    #第五步，利用freq_analysis_fasta_ID统计所有的fasta文件的频率（这一次VH直接用之前分析的去除头部冗余序列的fasta文件，所以不用重新分析）
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\conbined_VH_VL"
    #freq_analysis_fasta_ID(all_path)

    
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\conbined_VH_VL\edit_id"
    #split_char="-"
    #pos_start=3
    #pos_end=4
    #rm_file="T"
    #rm_repeat_seq(all_path,rm_file)#移除重复序列这样便于做进化树
    #split_id(all_path,split_char,pos_start,pos_end,rm_file)
    #phylo_draw(all_path)
    #第六步，根据已有的已经含有各个CDR的fasta文件，进行CDR3提取和频率分析
    









    ######第四次运行2018.11.15
     #第一步利用FGS模块将核苷酸序列合并。并且拆分为VH VL 两个文件，同时记录下来含有终止密码子的序列
    #提前定义好特异性的参数
    #file_path=r"E:\Program\project\yanli\FGS_2018.11.15\20181012-CD-Smix1234"
    #F=["CAGTTTTGGCACAGGCGGCC","CAGAATTGGCCCAGGCGGCC","CAGTTTTGGCACACGCGGCC"]
    #R=["GGGCCGTCGGTGGTC","GGGCCGTCGTGGTCG","GGCCGTCGTGGTCGA","GGCCGTCGGTGGTCG","GGGCCGTTGGTGGTCGACCT","GGCCCGTCGGTGGTCGACCT"]
    #split_all=["RAGGGGSGGGGSGGLGGAAV","RAGGGGSGGGGSGGLGGAAA","RAVGGGSGGGGSGGLGGAAA"]#这里需要在循环过程中不断地添加新的linker序列


    #正式进行分析
    #FGS_analysis(file_path,F,R)
    #record_abnormal(file_path,F,R)
    #file_path=r"E:\Program\project\yanli\FGS_2018.11.8\20181029seq"
    #FGS_analysis(file_path,F,R)
    #record_abnormal(file_path,F,R)

    #第二步利用得到的VH和VL序列在linux系统中利用hmmscan进行分析


    #第三步利用match_fasta_txt匹配路径中的fasta和txt的匹配，并且利用匹配到的结果进行分析，得到CDR区域的序列
    ######################VH
    #path=r"E:\Program\project\yanli\FGS_2018.11.15\hmm_out\longest_CDR\VH"
    #scheme="Longest_CDR_H"
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)
    #######################VL
    #path=r"E:\Program\project\yanli\FGS_2018.11.15\hmm_out\longest_CDR\VL"
    #scheme="Longest_CDR_L"
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)
    #为了下一步的进行需要新建文件夹并且将所有的VH和VL的CDR合并
    #os.chdir(path)
    #path_to_mk=os.path.split(os.sep)
    #os.mkdir("")
    #第四步利用merge_VH_VL将VH和VL的序列的CDR区域合并
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.15\hmm_out\longest_CDR\combined_VH_VL"
    #merge_VH_VL(all_path)

    #第五步，利用freq_analysis_fasta_ID统计所有的fasta文件的频率（这一次VH直接用之前分析的去除头部冗余序列的fasta文件，所以不用重新分析）
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.15\hmm_out\longest_CDR\combined_VH_VL"
    #freq_analysis_fasta_ID(all_path)

    
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.15\hmm_out\longest_CDR\combined_VH_VL\edit_id"
    #split_char="T7"
    #pos_start=-1
    #pos_end=None
    #rm_file="T"
    #rm_repeat_seq(all_path,rm_file)#移除重复序列这样便于做进化树
    #split_id(all_path,split_char,pos_start,pos_end,rm_file)
    #phylo_draw(all_path)
    #第六步，根据已有的已经含有各个CDR的fasta文件，进行CDR3提取和频率分析
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out\Longest_CDR\combined_VH_VL\CDR3_analysis"
    #rm_file="T"
    #get_cdr3(all_path,rm_file)
    #freq_analysis_fasta_ID(all_path)
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out\Longest_CDR\combined_VH_VL\CDR3_analysis"
    #第七步，合并所有的重链序列，然后取出CDR3,在进行频率分析
    all_path=r"E:\Program\project\yanli\combined_analysis(2018.11.8+2018.11.15)"
    combined_fasta(all_path)
    rm_file="T"
    get_cdr3(all_path,rm_file)
    freq_analysis_num_id(all_path)
    #第八步，对文件夹中所有的氨基酸序列进行比对分析
    clustal_all_folder(all_path)
    
    
    





























    
    #第二次运行，这次主要是统计CDRH最长的时候的频率
    #path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out\Longest_CDR\VH"
    #scheme="Longest_CDR_H"
   
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)

    #第二次运行，这次主要是统计CDRL最长的时候的频率
    #path="E:\Program\project\yanli\FGS_2018.11.8\hmm_out\Longest_CDR\VL"
    #scheme="Longest_CDR_L"
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)

    #第二次运行，这一次主要是合并VH VL的输出文件，然后把CDR进行合并（这一步在合并VH和VL的文件的时候可能需要手动操作）
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\conbined_VH_VL"
    #merge_VH_VL(all_path)


    #第二次运行，统计所有的fasta文件的频率（这一次VH直接用之前分析的去除头部冗余序列的fasta文件，所以不用重新分析）
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\conbined_VH_VL"
    #freq_analysis_fasta(all_path)
    #path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out_2\longest_CDR\conbined_VH_VL\edit_id_freq"
    
    #split_id(path)
    #freq_analysis_fasta_ID(path)



    ################
    ###    function1  #### 从fasta文件做进化树，所用到的fasta文件适合于删除重复序列以后的fasta文件的分析
    ###############
    #接下来进行进化树的左幅首先进行比对（调用align_own）
    #path="E:\Program\project\yanli\FGS_2018.11.8\hmm_out\phylo_test"
    #rm_file="T"
    #rm_repeat_seq(path,rm_file)#首先移除fasta文件中重复的序列
    #sep_char="-"#如果序列的ID较长，无法较好的放置在进化树中
    #pos_start=3#开始和结束位置确定了id的分割位置
    #pos_end=4
    #rm_file="T"
    #split_id(path,sep_char,pos_start,pos_end,rm_file)
    #phylo_draw(path)






    
    """
    #第一次运行，这时候还未把这一部分写作模块，因此比较多
    all_path=["E:\Program\project\yanli\FGS_2018.11.8\hmm_out"]
    for file_path in all_path: #对目标路径进行遍历
        #取出所有文件夹中的txt文件和fasta文件，并且成对放在字典中
        os.chdir(file_path)
        files_all=os.listdir(file_path)
        dict_txt_fasta={}
        for file_n in files_all:
            file_name,file_type=os.path.splitext(file_n)
            fasta_name=file_name+"."+"fasta"
            print fasta_name
            print file_name
            print file_type
            if file_type==".txt" and files_all.count(fasta_name)==1:
                dict_txt_fasta[file_n]=fasta_name#将配对成功的文件储存在字典中，方便以后调用
            print dict_txt_fasta
        for txt_file,fasta_file in dict_txt_fasta.items():
            file_name=os.path.splitext(txt_file)[0]
            file_cdr_seq=open("%s_CDR_sequences.fasta"%file_name,"w")
            hmmer_output=txt_file
            fasta_file=fasta_file
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
            print domain_sequences
            for query_id,domain_seq in domain_sequences.items():
                CDR1_seq=domain_seq["CDR1"]
                CDR2_seq=domain_seq["CDR2"]
                CDR3_seq=domain_seq["CDR3"]
                file_cdr_seq.write(">"+query_id+"\n")
                file_cdr_seq.write(CDR1_seq+"\t"+CDR2_seq+"\t"+CDR3_seq+"\n")
            file_cdr_seq.close()
            
    all_path="E:\Program\project\yanli\FGS_2018.11.8\hmm_out"
    merge_VH_VL(all_path)

    #这里分析fasta文件中序列的频率
    all_path="E:\Program\project\yanli\FGS_2018.11.8\hmm_out"
    freq_analysis_fasta(all_path)
    #清除重链序列首端的杂序列
    all_path="E:\Program\project\yanli\FGS_2018.11.8\hmm_out"
    head_str="AAA"
    remove_header_vh(all_path,head_str)


#第二次分析，将所有的CDR区域各种编号中最长的序列取出来，并且进行频率计算
#首先第一步进行CDR拆分
path=r"E:\Program\project\yanli\FGS_2018.11.8\hmm_out\Longest_CDR\VH"
"""

