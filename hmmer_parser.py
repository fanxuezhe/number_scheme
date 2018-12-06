#coding=utf-8
from __future__ import division
from Bio.SearchIO.HmmerIO import Hmmer3TextParser
import os
from Bio import SeqIO
from os import path,access,R_OK


##########################
#####FGS_ANALYSIS#########
##########################

#定义一个function,该function的功能是通过给定路径和引物（或者是需要的序列的前面的序列）以及VH，VL之间的linker序列，提取一代测序结果中的核苷酸序列，氨基酸序列，并且根据将其分成VH,VL两个序列，并且分别保存在不同的文件夹中
#输入的参数：路径，引物（正向和反向），链的类型(VH/VL),linker序列,rev_c 代表序列是否是反向测序
#输出的结果：保存在DNA_seq，AA_seq,VL_seq,VH_seq中的分别是带有引物的核苷酸序列，氨基酸序列的全长，VL和VH序列
def FGS_analysis(file_path,primer_f,primer_r,chain_type=None,rev_c=None,split_all=None):
    import os,re
    from Bio.Seq import Seq
    os.chdir(file_path)
    path_name=file_path.split(os.sep)[-1]
    start_pattern=re.compile("[EV][V][Q][LP]")#定义重链开始的正则表达式，定义这个是因为有部分的linker序列很短，但是数量也比较多，因此通过限制linker和气候的重链的开始序列来进行分析
    short_seq=open("sequences_too_short.fasta","w")
    DNA_seq=open("%s_sequences_DNA.fasta"%path_name,"w")
    AA_seq=open("%s_sequences_AA.fasta"%path_name,"w")
    if file_path.lower().find("fab")!=-1:
        chain_type=None
    if chain_type:
        if chain_type=="VH":
            VH_seq=open("%s_sequences_VH.fasta"%path_name,"w")
        elif chain_type=="VL":
            VL_seq=open("%s_sequences_VL.fasta"%path_name,"w")
    elif not chain_type:
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
            
            if rev_c:#判断是否是反向测序，假如是反向测序则将其调整成正常状态
                sequence_all=Seq(sequence_all).reverse_complement()
            elif not rev_c:
                pass
            for f_n in primer_f:
                if sequence_all.find(f_n)!=-1:
                    sequence_all_1=sequence_all.split(f_n)[1]
                    #print sequence_all_1
                    for r_n in primer_r:
                        print r_n
                        if sequence_all_1.find(r_n)!=-1:
                            sequence_all_2=sequence_all_1.split(r_n)[0]
                            aa_seq=str(Seq(sequence_all_1.split(r_n)[0]).translate())
                        
                            break
            #添加一种新的情况，假如序列为反向并且数量很少的时候可以直接使用接下来的过程进行处理
            if not aa_seq:
                seqence_all=Seq(sequence_all).reverse_complement()
                for f_n in primer_f:
                    if sequence_all.find(f_n)!=-1:
                        sequence_all_1=sequence_all.split(f_n)[1]
                        #print sequence_all_1
                        for r_n in primer_r:
                            print r_n
                            if sequence_all_1.find(r_n)!=-1:
                                sequence_all_2=sequence_all_1.split(r_n)[0]
                                aa_seq=str(Seq(sequence_all_1.split(r_n)[0]).translate())
                
            assert aa_seq,"there are primers that need to be added"
            """if file_path.lower().find("fab")!=-1:#增加一个判断，假如是fab文库那么直接将其中的氨基酸取出来并且按照相印的测序引物对应着将氨基酸序列写入相应的fasta文件中
                print aa_seq
                if file_n.upper().find("T7")!=-1:#这是VL的序列
                    name_id=file_n.upper().split("T7")[0]
                    DNA_seq.write(">"+"VL"+name_id+"\n")
                    DNA_seq.write(sequence_all_2+"\n")
                    VL_seq.write(">"+name_id+"\n")
                    VL_seq.write(aa_seq+"\n")
                    vl_seq=aa_seq
                if file_n.upper().find("T2A")!=-1:#这是VL的序列
                    name_id=file_n.upper().split("T2A")[0]
                    DNA_seq.write(">"+"VH"+name_id+"\n")
                    DNA_seq.write(sequence_all_2+"\n")
                    VH_seq.write(">"+name_id+"\n")
                    VH_seq.write(aa_seq+"\n")
                    vh_seq=aa_seq
            """
            if split_all or split_all_short:
                if aa_seq.find("*")==-1:
                #print sequence_all
                    print aa_seq
                    print file_name
                    #将DNA序列写到fasta文件中
                    DNA_seq.write(">"+str(file_name)+"\n")
                    DNA_seq.write(sequence_all_2+"\n")
                    #将氨基酸序列写到fasta文件中
                    AA_seq.write(">"+str(file_name)+"\n")
                    AA_seq.write(aa_seq+"\n")
                    #将氨基酸序列分为两个部分，VL和VH
                    for sp_n in split_all:
                        if aa_seq.find(sp_n)!=-1:
                            vl_seq,vh_seq=aa_seq.split(sp_n)
                            name_id="".join(file_name.split(" "))#由于当id中存在空格的时候，hmmscan处理数据时会将id自动断开，取第一部分作为id,因此去除掉
                            print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX%d"%i
                            print vl_seq
                            print vh_seq
                            if len(vh_seq)>85 and len(vl_seq)>80:
                                VL_seq.write(">"+name_id+"\n")
                                VL_seq.write(vl_seq+"\n")
                                VH_seq.write(">"+name_id+"\n")
                                VH_seq.write(vh_seq+"\n")
                            elif len(vh_seq)<=85 or len(vl_seq)<=80:
                                short_seq.write(">"+name_id+"\n")
                                short_seq.write(vl_seq+"\t"+vh_seq)
                            break
                        elif aa_seq.find(sp_n)==-1:
                            vh_seq=None
                            vl_seq=None
                    if not vl_seq and not vh_seq:#判断短的linker进行拆分的，短的linker后边需要保证直接是重链的序列
                        for sp_n in split_all_short:
                            if aa_seq.find(sp_n)!=-1:
                                vl_seq,vh_seq=aa_seq.split(sp_n)
                                name_id="".join(file_name.split(" "))#由于当id中存在空格的时候，hmmscan处理数据时会将id自动断开，取第一部分作为id,因此去除掉
                                print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX%d"%i
                                print vl_seq
                                print vh_seq
                                if re.match(start_pattern,vh_seq[:4]) and len(vl_seq)>80 and len(vh_seq)>85:
                                    VL_seq.write(">"+name_id+"\n")
                                    VL_seq.write(vl_seq+"\n")
                                    VH_seq.write(">"+name_id+"\n")
                                    VH_seq.write(vh_seq+"\n")
                                    break
                                elif re.match(start_pattern,vh_seq[:4]) and len(vl_seq)<=80 and len(vh_seq)<=85:
                                    short_seq.write(">"+name_id+"\n")
                                    short_seq.write(vl_seq+"\t"+vl_seq+"\n")
                                    break
                                elif not re.match(start_pattern,vh_seq[:4]):
                                    vl_seq=None
                                    vh_seq=None
                    if not vl_seq and not vh_seq:#假如短的linker依然没有找到那么判断它是不是fab库，假如是fab则直接把它写入到VH文件或者VL文件，假如都不是，那么vh_seq和vl_seq都变成None
                        
                        if file_n.upper().find("T2A")!=-1 and len(aa_seq)<150:#这时候有可能出现文件夹没有用fab命名，但是文件夹中有一部分fab的数据
                            name_id=file_n.upper().split("T2A")[0]
                            VH_seq.write(">"+name_id+"\n")
                            VH_seq.write(aa_seq+"\n")
                            vh_seq=aa_seq#这里为vl_seq赋值，这样能够避免后续vl_seq和vh_seq都变成none导致assert语句报错
                        elif file_n.upper().find("T2A")!=-1 and len(aa_seq)<85:#假如这里的序列太短那么就记录下来
                            name_id=file_n.upper().split("T2A")[0]
                            short_seq.write(">"+name_id+"\n")
                            short_seq.write(aa_seq+"\t")
                        elif file_n.upper().find("T7")!=-1 and len(aa_seq)<150:#这时候有可能出现文件夹没有用fab命名，但是文件夹中有一部分fab的数据,小于150保证是fab的文件
                            name_id=file_n.upper().split("T7")[0]
                            VL_seq.write(">"+name_id+"\n")
                            VL_seq.write(aa_seq+"\n")
                            vl_seq=aa_seq#这里为vh_seq赋值，这样能够避免后续vl_seq和vh_seq都变成none导致assert语句报错
                        elif file_n.upper().find("T7")!=-1 and len(aa_seq)<80:#假如这里的序列太短那么就记录下来
                            name_id=file_n.upper().split("T7")[0]
                            short_seq.write(">"+name_id+"\n")
                            short_seq.write(aa_seq+"\t")
                        elif  len(aa_seq)>=150:
                            vh_seq=None
                            vl_seq=None
                    print vl_seq
                    print vh_seq
                    assert vl_seq or vh_seq,"there is linker need to add"
                    
                    #这里利用i来记录这是第几条序列，同时也用来表示fasta文件中序列ID的唯一性
                    i+=1
                    print i
                else:
                    print "there is sequence contain the *"
            elif not split_all:
                DNA_seq.write(">"+str(file_name)+"\n")
                DNA_seq.write(sequence_all_2+"\n")
                #将氨基酸序列写到fasta文件中
                AA_seq.write(">"+str(file_name)+"\n")
                AA_seq.write(aa_seq+"\n")
                    #将氨基酸序列分为两个部分，VL和VH
                if chain_type=="VH":
                    vh_seq=aa_seq
                    name_id="".join(file_name.split(" "))#由于当id中存在空格的时候，hmmscan处理数据时会将文件名中的空格去除掉作为id
                    print vh_seq
                    VH_seq.write(">"+name_id+"\n")
                    VH_seq.write(vh_seq+"\n")
                elif chain_type=="VL":
                    vl_seq=aa_seq
                    name_id="".join(file_name.split(" "))#由于当id中存在空格的时候，hmmscan处理数据时会将文件名中的空格去除掉作为id
                    print vl_seq
                    VL_seq.write(">"+name_id+str(i)+"\n")
                    VL_seq.write(vl_seq+"\n")

        
    DNA_seq.close()           
    AA_seq.close()
    short_seq.close()
    if chain_type=="VH":
        VH_seq.close()
    if chain_type=="VL":
        VL_seq.close()








    
#定义一个新的功能，这个程序能够记录所有含有终止密码子的氨基酸序列，这些氨基酸序列可能是因为测序的原因导致最后的终止密码子的出现
#输入的参数：路径，引物（正向和反向）
def record_abnormal(path,primer_f,primer_r):
    import os
    from Bio.Seq import Seq
    os.chdir(path)
    files_all=os.listdir(path)
    file_abnormal=open("sequence_with_stop_codon.fasta","w")
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
            if aa_seq.find("*")!=-1:#debug发现在分析过程中发现有部分序列只有头和尾但是也有上下游引物但是没有中间的序列
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
            CDR1_seq=domain_seq[0]["CDR1"]
            CDR2_seq=domain_seq[0]["CDR2"]
            CDR3_seq=domain_seq[0]["CDR3"]
            species_n=domain_seq[1]
            ctype=domain_seq[2]
            evalue=domain_seq[3]
            file_cdr_seq.write(">"+query_id+"_"+species_n+"_"+ctype+"_"+str(evalue)+"\n")
            file_cdr_seq.write(CDR1_seq+"\t"+CDR2_seq+"\t"+CDR3_seq+"\n")#将CDR序列写入文件中
        file_cdr_seq.close()

#定义一个function，该function的功能是通过调用get_pos_domain等返回一个{id:{fr1:aaa,fr2:aaa....}}这样的格式的数据
        
def return_seq_domain(file_path,dict_txt_fasta,scheme="IMGT"):
    os.chdir(file_path)
    for txt_file,fasta_file in dict_txt_fasta.items():#对配对的fasta和txt文件进行遍历
        file_name=os.path.splitext(txt_file)[0]#获得txt文件的文件名
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
        return domain_sequences

            
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
                species, ctype = hsp.hit_id.split("_")
                evalue=hsp.evalue
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
        anno_domain.append((species,ctype,evalue,query_id,hmm_start,hmm_end,query_start,query_end,anno_domain_part))#最后一项为插入缺失位置和长度的记录，假如结构域为两个那么最后一项为多个列表的元组
        #print anno_domain
        print "this is %s"%anno_domain
    anno_seq=anno_domain
    return anno_seq

#function4:结合hmmer柄文件分析的结果,得出不同结构域的起始和结束的位置

def get_pos_domain(anno_seq,seq_dict,scheme="IMGT"):#目前只考虑一个domain的情况
    schemes_limit={"IMGT":[0,26,38,55,65,104,117],
            "AHO":[0,25,42,58,77,92,138],
            "Kabat_H":[0,31,40,54,74,106,117],
            "Kabat_K":[0,23,40,55,69,104,117],
            "Kabat_L":[0,23,40,55,69,104,117],
            "Chothia_H":[0,26,35,57,63,107,116],
            "Chothia_K":[0,25,32,55,58,106,116],
            "Chothia_L":[0,24,38,55,58,106,116],
            "LONGEST_CDR_H":[0,26,40,54,74,104,117],
            "LONGEST_CDR_L":[0,23,40,55,69,104,117] }
    scheme_limit=schemes_limit[scheme.upper()]
    region_anti=["FR1","CDR1","FR2","CDR2","FR3","CDR3"]
    region_seq={}#用来记录每个序列的各个结构域的序列
    for seq_anno in anno_seq:
        species_n=seq_anno[0]
        print species_n
        ctype=seq_anno[1]
        print ctype
        evalue=seq_anno[2]
        print evalue
        hmm_start=seq_anno[4]#hmm_start,hmm_end,query_start,query_end在目前仅考虑只有一个结构域的情况
        hmm_end=seq_anno[5]
        query_start=seq_anno[6]
        query_end=seq_anno[7]
        if len(seq_anno[-1])>1:#取出含有插入和缺失位置信息的列表，列表内容为'M03186:269:000000000-BK7T5:1:1101:17338:2018:CCTCCTGA+TCTTTCCC', 0, 112, 0, 121, [([(98, 9)], [])]的最后一项
            domain_details=seq_anno[-1][0]
        elif len(seq_anno[-1])==1:
            domain_details=seq_anno[-1]
        domain_num=1#目前仅遍历第一个结构域
        insert_num=len(domain_details[0][0])
        del_num=len(domain_details[0][1])#domain[0代表取出第一个domain]
        #domain_num=len(domain_details)
        len_to_add=min(hmm_start,query_start)#确定最小的开始位置，以免后续获取序列的时候,索引出现问题
        if hmm_start<10:
            seq_start=query_start-hmm_start#找到序列开始位置的索引，这也是FR1开始的位置
        elif hmm_start>10:
            pass#暂时还不知道如何处理
        
        query_id=seq_anno[3]#这里来处理query的id
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
        region_seq[query_id]=(domain_seq,species_n,ctype,evalue)
    return region_seq
def merge_VH_VL(file_path):#这部分功能主要是将VH和VL的CDR区进行合并随后在进行一个频率分析,species为将要排除的物种
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
        #assert len(vl_all)==len(vh_all),"attention the number of CDR of VL is different from the numner of CDR of VH"
        dict_vl={}
        dict_vh={}
        #这个循环的目的是把VH和VL的id和序列都储存在字典中，这样就可以方便后续用ID进行匹配序列
        for i in range(0,len(vl_all),2):
            print vl_all[i]
            print vl_all[i+1]
            if vl_all[i][0]==">" and vl_all[i+1][0]!=">":
                #species=vl_all[i].split("_")[-3]
                dict_vl[vl_all[i]]=vl_all[i+1]
        for i in range(0,len(vh_all),2):
            print vh_all[i]
            print vh_all[i+1]
            if vh_all[i][0]==">" and vh_all[i+1][0]!=">":
                dict_vh[vh_all[i]]=vh_all[i+1]
            else:
                print "the fasta file is not standard the sequence not follow the id"
        #接下来要合并两个文件,文件合并时分为两种情况第一种都匹配（scFv），第二种两者有部分不匹配(Fab),这种情况下应该取出相同的放在一个文件中，不相同的分别放在不同的文件中
        file_name=os.path.splitext(vl_cdr)[0]
        combine_cdr=open("%s_combined_VH_VL_CDR.fasta"%file_name,"w")
        vl_unique=open("%s_VL_unique.fasta"%file_name,"w")
        vh_unique=open("%s_VH_unique.fasta"%file_name,"w")
        common_id=[]
        print dict_vl
        print dict_vh
        for vl_id,vl_seq in dict_vl.items():#对VL中的序列进行遍历，并且找出两者共有的序列，并且取出ID
            if dict_vh.has_key(vl_id):#首先检查VH中是否有对应的ID,假如不存在则写入VL_unique相应的文件中
                combined_seq=vl_seq.strip("\n")+"\t"+dict_vh[vl_id]
                combine_cdr.write(vl_id)
                combine_cdr.write(combined_seq)
                common_id.append(vl_id)
            else:
                vl_unique.write(vl_id)
                vl_unique.write(vl_seq)
        for vh_id,vh_seq in dict_vh.items():
            if common_id.count(vh_id)==0:#当VH的字典中没有相应的id的时候，则记录下来
                vh_unique.write(vh_id)
                vh_unique.write(vh_seq)
                
        combine_cdr.close()
        vh_unique.close()
        vl_unique.close()
        

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
def freq_analysis_fasta_ID(path,len_out):#这里是指将来在文件中的第一行注释将会有几列，当len_out=1的时候，默认认为是CDRH3,当len_out=6的时候默认为L+H两部分组成
    import os,shutil
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta" and file_name.find("freq")==-1:
            all_seq=open(file_n,"r").readlines()
            total_num=(len(all_seq))/2
            unique_seq=[]
            dict_unique_seq={}#利用一个字典对唯一序列进行记录这样方便后续对ID的书写
            freq_seq={}
            header=file_n.split(".")[0]
            file_freq=open("%s_freq.fasta"%header,"w")
            if len_out==1:
                file_freq.write("出现次数"+"\t"+"H-CDR3"+"\t"+"IDs"+"\n")
            elif len_out==6:
                file_freq.write("出现次数"+"\t"+"L-CDR1"+"\t"+"L-CDR2"+"\t"+"L-CDR3"+"\t"+"H-CDR1"+"\t"+"H-CDR2"+"\t"+"H-CDR3"+"\t"+"IDs"+"\n")
            elif len_out==3:
                file_freq.write("出现次数"+"\t"+"H-CDR1"+"\t"+"H-CDR2"+"\t"+"H-CDR3"+"\t"+"IDs"+"\n")
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
                ids=[x.strip("\n").strip(">")  for x in freq_n[1]]
                file_freq.write(str(freq)+"\t"+str(seq_n).strip("\n")+"\t"+"  ".join(ids)+"\n")#2018.11.23在输出的文件中删掉“>”
            file_freq.write(str(len(unique_seq))+"\t"+str(len(unique_seq))+"/"+str(total_num)+"="+str(len(unique_seq)/total_num))#输出相应的比率
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

#这里定义一个function用来从A文件中提取出B文件对应的核苷酸序列，并且切除头尾的序列，

from Bio.Seq import Seq
def get_right_DNA(path,file_1,file_2,primer_pair1,primer_pair2,len_del):#file_1是氨基酸文件，file_2是DNA序列文件，其他参数在下方的解释中
    os.chdir(path)#转到目标文件夹中
    file_1_content=open(file_1,"r").readlines()#读取文件的内容
    file_2_content=open(file_2,"r").readlines()
    file_1_dict={}#定义空字典将来用于储存id和序列，并且一一对应
    file_2_dict={}
    file_name=file_2.split(".")[0]#找到除了后缀以外的文件名称，便于未来命名
    file_dna=open("%s_DNA_seq.fasta"%file_name,"w")#未来的输出文件的名称
    for i in xrange(len(file_1_content)):
        if file_1_content[i].startswith(">"):#确认该位置是不是id
            file_1_dict[file_1_content[i]]=file_1_content[i+1]#将序列储存在字典中，键为id，值为序列
    for i in xrange(len(file_2_content)):
        if file_2_content[i].startswith(">"):
            file_2_dict[file_2_content[i]]=file_2_content[i+1]
    dna_seq_formal_1=""
    dna_seq_formal=""
    for id_n,seq_n in file_1_dict.items():
        dna_seq=file_2_dict[id_n]
        for f,r in primer_pair1.items():
            if dna_seq.find(f)!=-1 and dna_seq.find(r)!=-1:
                dna_seq_formal="".join(("".join(dna_seq.split(f)[1:])).split(r)[:-1])[:-len_del]
        for f,r in primer_pair2.items():
            if dna_seq.find(f)!=-1 and dna_seq.find(r)!=-1:
                dna_seq_formal="".join(("".join(dna_seq.split(f)[1:]).split(r)[:-1]))[len_del:]
                dna_seq_formal=str(Seq(dna_seq_formal.strip("\n")).reverse_complement())
        if not dna_seq_formal.endswith("\n"):
            dna_seq_formal=dna_seq_formal+"\n"
        file_dna.write(id_n)
        file_dna.write(dna_seq_formal)
        
    file_dna.close()


"""
这里需要提供的参数的例子有：
#第一个文件夹        
path="E:\Program\project\yanli\NGS_1st\macaca_fascicularis\VH_100_VK95_new\Genesp_VH"
file_1="Genesp-VH_R1_out_merged_new_AA.fasta"#含有正常序列的氨基酸序列文件
file_2="Genesp-VH_R1_out_merged.fasta"#含有原始的核苷酸序列的fasta文件
primer_pair1={"AAGCGGCCGCT":"GTCTTCCC"}#第一对引物，这里以字典的形式存在，它的内容为正向序列的头尾引物
primer_pair2={"GGGAAGAC":"AGCGGCCGC"}#第二对引物，这里以字典的形式存在，他的内容为反向序列的头尾引物
len_del=21#这里用来确定尾部序列需要去除的序列长度
get_right_DNA(path,file_1,file_2,primer_pair1,primer_pair2,len_del)

"""
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
#接下来定义一个新的function，该功能是利用已有的含有CDR1/2/3的fasta文件转化成为含有CDR3的序列,序列的ID不变,仅处理fasta文件
def get_cdr3(path,rm_file="F"):
    os.chdir(path)
    files_all=os.listdir(path)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta" and file_name.find("CDR3_only")==-1:#保证当重复计算时不会出现对CDR3_only文件进行分析的情况
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
#定义一个function，该function的功能是合并多个fasta文件，即目标文件夹下所有的fasta文件
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
#定义一个function，该function的功能是根据返回的FR和CDR序列就进行编号
"""
def number_aho(fr_cdr):
    fr1_seq=fr_cdr["FR1"]
    cdr1_seq=fr_cdr["CDR1"]
    fr2_seq=fr_cdr["FR2"]
    cdr2_seq=fr_cdr["CDR2"]
    fr3_seq=fr_cdr["FR3"]
    cdr3_seq=fr_cdr["CDR3"]
    #首先检测FR1的前十位有没有插入或者缺失，前十位的插入缺失将会安排在第8的位置
    len_fr1=len(fr1_seq)
    #if len_fr1>
    fr1_num=
    cdr_num=
    fr2_num=
    cdr2_num=
    fr3_num=
    cdr3_num=
def imgt_cdr(cdr_seq):



def(start_pos,end_pos,mid
    
"""
#定义一个function，该function功能为将列表中所有的元素都变为列表
def transfer_all_str(list_n):
    print list_n
    list_m=[str(i) for i in list_n]
    print list_m
    return list_m
def number_imgt(fr_cdr):
    all_cdr_numbered=[]
    char="\t"
    for seq_id,sequences in fr_cdr.items():
        fr1_seq=sequences["FR1"]
        cdr1_seq=sequences["CDR1"]
        fr2_seq=sequences["FR2"]
        cdr2_seq=sequences["CDR2"]
        fr3_seq=sequences["FR3"]
        cdr3_seq=sequences["CDR3"]
        alphabet=["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
        #首先检测FR1的前十位有没有插入或者缺失，前十位的插入缺失将会安排在第8的位置
        #len_fr1=len(fr1_seq)
        #if len_fr1>
        
        #首先判断CDR1的长度，假如超过最大长度（12），那么在两侧分布，假如不足12则在两侧分别减少    
        cdr1_len=len(cdr1_seq)
        if cdr1_len<=12:
            _len=12-cdr1_len
            char_add="-"*(_len)
            if cdr1_len%2:
                left_len=cdr1_len/2+1
                right_len=cdr1_len/2
            elif not cdr1_len%2:
                left_len=cdr1_len/2+1
                right_len=cdr1_len/2
            left_seq=cdr1_seq[:(left_len-1)]
            right_seq=cdr1_seq[-(right_len):]
            cdr1_num=char.join(transfer_all_str(range(27,27+left_len+1)))+char_add+char.join(transfer_all_str(range(27+(12-right_len),39)))#13是指12加一
            cdr1_seq_all=left_seq+char_add+right_seq
        if cdr1_len>12:
            if cdr1_len%2:
                left_len=cdr1_len/2+1
                right_len=cdr1_len/2
            elif not cdr1_len%2:
                left_len=cdr1_len/2
                right_len=cdr1_len/2
            left_seq=cdr1_seq[:(left_len-1)]
            right_seq=cdr1_seq[-(right_len):]
            cdr1_num=char.join(transfer_all_str(range(27,32)))+char.join(zip((["32"]*(cdr1_len-12)),alphabet[:cdr1_len-12-1]))+char.join(transfer_all_str(range(33,38+1)))
            cdr1_seq_all=left_seq+right_seq
        #对CDR2进行编号
        cdr2_len=len(cdr2_seq)
        if cdr2_len<=10:
            _len=10-cdr2_len
            char_add="-"*(_len)
            if cdr2_len%2:
                left_len=cdr2_len/2+1
                right_len=cdr2_len/2
            elif not cdr2_len%2:
                left_len=cdr2_len/2
                right_len=cdr2_len/2
            left_seq=cdr2_seq[:(left_len-1)]
            right_seq=cdr2_seq[-(right_len):]
            cdr2_num=char.join(transfer_all_str(range(56,56+left_len+1)))+char_add+char.join(transfer_all_str(range(56+(10-right_len),66)))#13是指12加一
            cdr2_seq_all=left_seq+char_add+right_seq
        if cdr2_len>10:
            if cdr2_len%2:
                left_len=cdr2_len/2+1
                right_len=cdr2_len/2
            elif not cdr2_len%2:
                left_len=cdr2_len/2+1
                right_len=cdr2_len/2
            left_seq=cdr2_seq[:(left_len-1)]
            right_seq=cdr2_seq[-(right_len):]
            cdr2_num=char.join(transfer_all_str(range(56,60)))+char.join(zip((["60"]*(cdr2_len-10)),alphabet[:cdr2_len-10-1]))+char.join(transfer_all_str(range(61,65+1)))
            cdr2_seq_all=left_seq+right_seq
        #对CDR3进行编号
        cdr3_len=len(cdr3_seq)
        if cdr3_len<=13:
            _len=10-cdr1_len
            char_add="-"*(_len)
            if cdr3_len%2:
                left_len=cdr3_len/2+1
                right_len=cdr3_len/2
            elif not cdr3_len%2:
                left_len=cdr3_len/2
                right_len=cdr3_len/2
            left_seq=cdr3_seq[:(left_len-1)]
            right_seq=cdr3_seq[-(right_len):]
            cdr3_num=char.join(transfer_all_str(range(105,105+left_len+1)))+char_add+char.join(transfer_all_str(range(105+(13-right_len),118)))#13是指12加一
            cdr3_seq_all=left_seq+char_add+right_seq
        if cdr3_len>13:
            if cdr3_len%2:
                left_len=cdr3_len/2+1
                right_len=cdr3_len/2
            elif not cdr3_len%2:
                left_len=cdr3_len/2+1
                right_len=cdr3_len/2
            left_seq=cdr3_seq[:(left_len-1)]
            right_seq=cdr3_seq[-(right_len):]
            if (cdr3_len-13)%2:#检查CDR3比13长出来多少，假如是个奇数则在N端多分配一个氨基酸，否则两侧均匀分配
                part_1=char.join(transfer_all_str(range(105,111+1)))#加一是因为range少1
                a="111"
                part_2=char.join([(a+"."+str(b)) for b in range(1,((cdr3_len-13)/2+2))])
                print part_2
                a="112"
                part_3=char.join([(a+"."+str(b)) for b in range(((cdr3_len-13)/2),0,-1)])
                print part_3
                part_4=char.join(transfer_all_str(range(112,(117+1))))
                cdr3_num=part_1+"\t"+part_2+"\t"+part_3+"\t"+part_4
            elif not (cdr3_len-13)%2:
                part_1=char.join(transfer_all_str(range(105,111+1)))
                a="111"
                part_2=char.join([(a+"."+str(b)) for b in range(1,((cdr3_len-13)/2)+1)])
                print part_2
                a="112"
                part_3=char.join([(a+"."+str(b)) for b in range(((cdr3_len-13)/2),0,-1)])
                print part_3
                part_4=char.join(transfer_all_str(range(112,(117+1))))
                cdr3_num=part_1+"\t"+part_2+"\t"+part_3+"\t"+part_4
            cdr3_seq_all=left_seq+right_seq
        seq_numbered={seq_id:{"CDR1":{cdr1_num:cdr1_seq_all},"CDR2":{cdr2_num:cdr2_seq_all},"CDR3":{cdr3_num:cdr3_seq_all}}}
        all_cdr_numbered.append(seq_numbered)
    return all_cdr_numbered 


"""def number_kabat_h(fr_cdr):
    fr1_seq=fr_cdr["FR1"]
    cdr1_seq=fr_cdr["CDR1"]
    fr2_seq=fr_cdr["FR2"]
    cdr2_seq=fr_cdr["CDR2"]
    fr3_seq=fr_cdr["FR3"]
    cdr3_seq=fr_cdr["CDR3"]
    #首先对CDR1进行编号
    alphabet=["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]
    char="-"
    len_cdr1=len(cdr1_seq)#cdr1插入规则：直接插在35后方，缺失规则从35-31逐步缺失
    if len_cdr1<=5:
        _len=5-len_cdr1
        char_add="-"*_len
        cdr1_seq_all=cdr1_seq+char_add
        cdr1_num=range(31,36)
    elif len_cdr1>5:
        cdr1_seq_all=cdr1_seq
        start_char="35"
        part_1=transfer_all_str(range(31,35+1))
        part_2=[(start_char+alphabet[i])for i in range(len(len_cdr1-5))]
        cdr1_num=part_1+"\t"+part_2
    len_cdr2=len(cdr2_seq)#cdr2插入规则，直接插入在52后方，缺失规则：先缺失52/51/50再缺失53-65
    if len_cdr2<=16:#当CDR2没有插入的时候,这时候编号不用变
        _len=16-len_cdr2
        to_add_char="-"*_len
        if _len<3:
            cdr2_num=char.join(transfer_all_str(range(50,65+1)))
            cdr2_seq_all=cdr2_seq[:3-_len]+to_add_char+cdr2_seq[3-_len：]
        elif _len>=3:
            cdr2_num=char.join(transfer_all_str(range(50,65+1)))
            cdr2_seq_all=to_add_char+cdr2_seq[_len:]
    elif len_cdr>16:
        _len=len_cdr2-16
        start_char="52"
        to_add_char=[(start_char+alphabet[i]) for i in range(_len)]
        cdr2_num=transfer_all_str(range(50,52+1))+to_add_char+transfer_all_str(range(53,65))
        cdr2_seq_all=cdr2_seq
            
        
   len_cdr3=len(cdr3_seq)#cdr2插入规则，直接插入在52后方，缺失规则：先缺失52/51/50再缺失53-65
    if len_cdr3<=8:#当CDR2没有插入的时候,这时候编号不用变
        _len=8-len_cdr3
        to_add_char="-"*_len
        if _len<6:
            cdr2_num=char.join(transfer_all_str(range(95,102+1)))
            cdr2_seq_all=cdr2_seq[:3-_len]+to_add_char+cdr2_seq[6-_len：]
        elif _len>=6:
            cdr2_num=char.join(transfer_all_str(range(50,65+1)))
            cdr2_seq_all=to_add_char+cdr2_seq[_len:]
    elif len_cdr>16:
        _len=len_cdr2-16
        start_char="52"
        to_add_char=[(start_char+alphabet[i]) for i in range(_len)]
        cdr2_num=transfer_all_str(range(50,52+1))+to_add_char+transfer_all_str(range(53,65))
        cdr2_seq_all=cdr2_seq
            
    

def number_kabat_l(fr_cdr):
    fr1_seq=fr_cdr["FR1"]
    cdr1_seq=fr_cdr["CDR1"]
    fr2_seq=fr_cdr["FR2"]
    cdr2_seq=fr_cdr["CDR2"]
    fr3_seq=fr_cdr["FR3"]
    cdr3_seq=fr_cdr["CDR3"]
    #首先检测FR1的前十位有没有插入或者缺失，前十位的插入缺失将会安排在第8的位置
    len_fr1=len(fr1_seq)
    if len_fr1>
    fr1_num=
    cdr_num=
    fr2_num=
    cdr2_num=
    fr3_num=
    cdr3_num=


def number_chothia_h(fr_cdr):
    fr1_seq=fr_cdr["FR1"]
    cdr1_seq=fr_cdr["CDR1"]
    fr2_seq=fr_cdr["FR2"]
    cdr2_seq=fr_cdr["CDR2"]
    fr3_seq=fr_cdr["FR3"]
    cdr3_seq=fr_cdr["CDR3"]
    #首先检测FR1的前十位有没有插入或者缺失，前十位的插入缺失将会安排在第8的位置
    len_fr1=len(fr1_seq)
    if len_fr1>
    fr1_num=
    cdr_num=
    fr2_num=
    cdr2_num=
    fr3_num=
    cdr3_num=



def number_chothia_l(fr_cdr):
    fr1_seq=fr_cdr["FR1"]
    cdr1_seq=fr_cdr["CDR1"]
    fr2_seq=fr_cdr["FR2"]
    cdr2_seq=fr_cdr["CDR2"]
    fr3_seq=fr_cdr["FR3"]
    cdr3_seq=fr_cdr["CDR3"]
    #首先检测FR1的前十位有没有插入或者缺失，前十位的插入缺失将会安排在第8的位置
    len_fr1=len(fr1_seq)
    if len_fr1>
    fr1_num=
    cdr_num=
    fr2_num=
    cdr2_num=
    fr3_num=
    cdr3_num=
            
        
"""            
#这里定义一个function，该function的功能是根据给出的字符串的类型将文件分为几个组,该程序只用来分类fasta文件
def divide_group(path,char):
    import shutil
    os.chdir(path)
    files_all=os.listdir(path)
    #首先根据给出的字符分别建立文件夹
    for char_n in char:
        os.mkdir(char_n)
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".seq":#判断文件类型是符合
            print "okkkkkkk"
            for char_n in char:
                print char_n
                if file_n.find(char_n)!=-1:
                    path_n=os.path.join(path,char_n)
                    shutil.copy(file_n,path_n)
    
#定义一个新的功能，该功能能够将含有CDR序列的文件中的特殊物种的序列提取出来，并且重新命名
def divide_species_cdr(path,species):
    os.chdir(path)
    files_all=os.listdir(path)
    all_species=["human","rat","rabbit"]
    del_num_record=open("sequence_number_removed.txt","w")
    for file_n in files_all:
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta" and file_name.find("CDR_sequences")!=-1:
            file_content=open(file_n,"r").readlines()
            file_name_to_use=file_name.split("CDR_sequences")[0]
            file_species=open("%s_%s_CDR_sequences.fasta"%(file_name_to_use,species),"w")
            file_without_species=open("%s_without_%s_CDR_sequences.fasta"%(file_name_to_use,species),"w")
            id_seq_dict={}
            id_seq_withou_species={}
            for i in xrange(len(file_content)):#首先对文件内容进行遍历，将含有特定物种名称的序列提取出来并且保存在字典中
                if file_content[i].startswith(">") and file_content[i].find(species)!=-1:
                    id_seq_dict[file_content[i]]=file_content[i+1]
                elif file_content[i].startswith(">") and file_content[i].find(species)==-1:
                    id_seq_withou_species[file_content[i]]=file_content[i+1]
            num_rm=len(id_seq_dict)
            for seq_id,seq_n in id_seq_dict.items():
                file_species.write(seq_id)
                file_species.write(seq_n)
            for seq_id,seq_n in id_seq_withou_species.items():
                file_without_species.write(seq_id)
                file_without_species.write(seq_n)
            del_num_record.write(file_n+str(num_rm)+"\n")
            file_species.close()
            file_without_species.close()
    del_num_record.close()

#定义一个功能,这个程序能够通过分析VH,VL中含有特定物种的id，并且将前面的id记录下来，并且返回这些ids,同时记录下不同文件夹的
def get_species_id(path,species):
    os.chdir(path)
    files_all=os.listdir(path)
    file_name_match={}
    for file_n in files_all:
        
        file_name,file_type=os.path.splitext(file_n)
        if file_type==".fasta":#判断文件是不是fasta文件
            if file_n.find("VL")!=-1 and file_n.startswith("C"):
                #print file_n
                file_vh=file_n.replace("VL","VH")
                #print file_vh
                #print files_all
                n=files_all.count(file_vh)
                #print n
                if files_all.count(file_vh)==1:
                    file_name_match[file_n]=file_vh
    print file_name_match
    species_id={}
    
    for file_vl,file_vh in file_name_match.items():#储存匹配文件的所有的fasta的id
        vl_name=os.path.splitext(file_vl)[0]
        file_combined_name=(vl_name+"_edited_id_combined_VH_VL_CDR.fasta")#这里定义合并后文件的名称方便后续的处理合并文件的时候的使用
        vl_content=open(file_vl,"r").readlines()
        vh_content=open(file_vh,"r").readlines()
        all_seq=vl_content+vh_content
        all_ids=[]
        for line_n in all_seq:
            if line_n.startswith(">") and line_n.find(species)!=-1:
                seq_id="".join(line_n.split(".")[:1])
                #print seq_id
                if all_ids.count(seq_id)==0:#去除重复
                    all_ids.append(seq_id)
        species_id[file_combined_name]=all_ids
        print file_combined_name
        print len(all_ids)
    return species_id

def rm_species_id(path,species,species_id):
    os.chdir(path)
    files_all=os.listdir(path)
    num_record=open("remove_num_record.txt","w")
    num_record.write("file_name"+"\t"+"removed_ids")
    for file_name,rm_ids in species_id.items():
        if files_all.count(file_name)==1:
            file_content=open(file_name,"r").readlines()
            file_name_vl=file_name.split("VL")[0]
            file_without_species=open("%s_without_%s_combined_VH_VL_CDR.fasta"%(file_name_vl,species),"w")
            
            x=0
            for i in xrange(len(file_content)-1):
                
                line_n=file_content[i]
                seq_n=file_content[i+1]
                #print file_content
                #print line_n
                if line_n.startswith(">") and rm_ids.count(line_n.strip("\n"))==0:
                    #print line_n
                    file_without_species.write(line_n)
                    file_without_species.write(seq_n)
                elif line_n.startswith(">") and rm_ids.count(line_n.strip("\n"))==1:
                    print line_n
                    x=x+1
            file_without_species.close()
            num_record.write(file_name+"\t"+str(x)+"\n")
    num_record.close()
                    
                        
#定义一个新建文件夹的功能，假如这个文件夹存在那么就删除掉并且新建新的，这样就可以反复的调试程序文件
def mkdir_folder(folder_path):#提供该需要新建文件夹的绝对路径
    import shutil,os
    if os.path.exists(folder_path):#判断两个文件夹是否存在，假如不存在那么就建立新的文件夹，假如存在那么就是删除后重新建立
        print "ooooook"
        shutil.rmtree(folder_path)
    os.mkdir(folder_path)

#定义一个功能，这个程序的功能是将给定的文件夹中的所有文件夹都列出来，接下来能够对其进行遍历
def get_folder(path):
    file_folder=os.listdir(path)
    all_folder=[]
    for file_n in file_folder:
        folder_n=(path+"\\"+file_n)
        if os.path.isdir(folder_n) and all_folder.count(folder_n)==0  and file_n!="result" and file_n!="fasta_file" and  file_n!="finished_folder":
            all_folder.append(folder_n)
    return all_folder

#定义一个功能，这个功能能够将现在文件夹中的正常VH和VL的序列拷贝到VH_VL氨基酸库的文件夹中,同时能够将含有终止密码子的fasta文件转移到对应的文件夹中
def cp_aa_fasta(folder_now,result_folder):#这里的target_folder应该是result_folder对应的文件夹
    os.chdir(folder_now)
    folder_name=os.path.split(folder_now)[-1]
    all_fasta=(result_folder+"\\fasta_file")
    abnormal_fasta=(result_folder+"\\"+folder_name+"\\abnormal_seq")
    target_folder=(result_folder+"\\"+folder_name+"\\AA_sequences")
    DNA_folder=(result_folder+"\\"+folder_name+"\\DNA_sequences")
    
    if os.path.exists("%s_sequences_VH.fasta"%folder_name) and os.path.getsize("%s_sequences_VH.fasta"%folder_name):
        shutil.copy("%s_sequences_VH.fasta"%folder_name,all_fasta)
        shutil.copy("%s_sequences_VH.fasta"%folder_name,target_folder)
    if os.path.exists("%s_sequences_VL.fasta"%folder_name) and os.path.getsize("%s_sequences_VL.fasta"%folder_name):
        shutil.copy("%s_sequences_VL.fasta"%folder_name,all_fasta)
        shutil.copy("%s_sequences_VL.fasta"%folder_name,target_folder)
    if os.path.exists("sequence_with_stop_codon.fasta") and os.path.getsize("sequence_with_stop_codon.fasta"):
        shutil.copy("sequence_with_stop_codon.fasta",abnormal_fasta)
    if os.path.exists("%s_sequences_DNA.fasta"%folder_name) and os.path.getsize("%s_sequences_DNA.fasta"%folder_name):
        shutil.copy("%s_sequences_DNA.fasta"%folder_name,DNA_folder)
    if os.path.exists("sequences_too_short.fasta") and os.path.getsize("sequences_too_short.fasta"):
        shutil.copy("sequences_too_short.fasta",abnormal_fasta)
  
#定义一个功能，这个功能能够通过给定的文件夹对文件中的测序文件进行分析，并且将已经完成的文件夹转移到finished_folder中去，这样能够避免反复的分析
def scfv_analysis(home_folder,all_folders):#home_folder就是所有一代测序文件夹所在的文件夹
    import shutil
    for i in xrange(len(all_folders)):
        folder_n=all_folders[i]
        folder_name=os.path.split(folder_n)[-1]#取出当前所在文件夹的文件夹名称，便于后续组织化的储存数据(比如异常序列的数据)
        os.chdir(folder_n)
        #随后进行翻译和异常序列的分析
        rev_c=None
        chain=None
        FGS_analysis(folder_n,F,R,chain,rev_c,split_all)
        record_abnormal(folder_n,F,R)
        cp_aa_fasta(folder_n,result_folder)#将正常的氨基酸序列拷贝到fasta文件库中去
        
        #然后将分析完的序列文件夹转移到finished_folder中去
        os.chdir(home_folder)
        folder_target=(home_folder+"\\"+"finished_folder"+"\\"+folder_name)
        shutil.move(folder_n,folder_target)
#定义一个功能，将文件夹A中的所有特定文件拷贝到B文件夹中去,type_char代表是一种字符类型，当其为VH时，在后边的文件夹中新建VH文件夹并且将含有VH的序列转移到这个VH文件夹中，这样便于后续的匹配和分析
def cp_all_file(folder_A,folder_B,type_char,mkdir="F"):
    import shutil,os
    files_all=os.listdir(folder_A)
    os.chdir(folder_A)
    folder_char=(folder_B+"\\"+type_char)
    if type_char and mkdir=="T" and not os.path.exsits(folder_char):
        mkdir_folder(folder_char)
    for file_n in files_all:
        if type_char and file_n.find(type_char)!=-1 and mkdir=="T" and os.path.isfile(file_n):#最后一个选项一定要保证要拷贝的文件不是文件夹
            shutil.copy(file_n,folder_char)
        elif type_char and file_n.find(type_char)!=-1 and mkdir=="F" and os.path.isfile(file_n):
            shutil.copy(file_n,folder_B)
#接下来定义一个功能，能够将一个文件从A文件夹拷贝到B文件夹（同时考虑到文件的名称的开始和结尾的位置是否符合特定的要求）
def cp_start_end(folder_A,folder_B,pos_char,pos="start"):
    import shutil,os
    files_all=os.listdir(folder_A)
    os.chdir(folder_A)
    for file_n in files_all:
        if pos=="start":
            if file_n.startswith(pos_char):
                shutil.copy(file_n,folder_B)
        elif pos=="end":
            if file_n.endswith(pos_char):
                shutil.copy(file_n,folder_B)

    
        
#为result文件夹进行组织,在result文件夹中为每一个储存测序文件的文件夹新建一个文件夹，并且每一个文件夹中都将包含四个文件夹（DNA_sequences,abnormal_seq,CDR_frequency,CDR3_frequency）
def organize_result(file_folder,folders):#这里的file_folder指的是各个含有测序文件的文件夹所在的文件夹
    finished_folder=(file_folder+"\\finished_folder")
    result_folder=(file_folder+"\\result")
    all_fasta=(result_folder+"\\fasta_file")
    mkdir_folder(result_folder)
    mkdir_folder(finished_folder)
    mkdir_folder(all_fasta)
    for folder_n in folders:
        folder_name=os.path.split(folder_n)[-1]
        folder_1=(result_folder+"\\"+folder_name)
        mkdir_folder(folder_1)
        for folder_m in ["AA_sequences","DNA_sequences","abnormal_seq","CDR_frequency","CDR3_frequency"]:#在这里定义每个文件夹中需要添加的文件夹名称
            folder_2=(folder_1+"\\"+folder_m)
            mkdir_folder(folder_2)
    
    

#定义一个新的功能，利用其记录一个文件中特定物种的抗体ID，并且记录其总的数目
def record_anti_id(path,species):
    os.chdir(path)
    file_species=open("sequeces_id_contain_%s"%species,"w")
    files_all=os.listdir(path)
    for file_n in files_all:
        if file_n.endswith(".fasta"):
            file_species.write(file_n+"\n")
            with open(file_n,"r") as all_seq:
                sequences=all_seq.readlines()
                i=0
                all_seq_id=[]
                for i in xrange(len(sequences)):
                    test_seq_id=sequences[i]
                    test_seq_seq
                    if test_seq_id.startswith(">") and test_seq_id.find(species)!=-1:
                        i=i+1
                        all_seq_id.append(test_seq_id)
                        all_seq_id.append(test_seq_seq)
                file_species.write("there is %d sequences containing the %s species"%(i,species))
                file_species.write(all_seq_id)
    file_species.close()
                        
#定义一个function将分析以后的结果拷贝到对应的文件夹中
def cp_file_to_result(file_repo):
    C_CDR3_folder=(file_repo+"\\all_match_file\\C_file\\C_CDR3\\formal_id")
    VH_CDR3_folder=(file_repo+"\\all_match_file\\VH\\VH_CDR3\\formal_id")
    VH_VL=(file_repo+"\\all_match_file\\VH_VL_file\\combined_VH_VL\\all_combined")
    C_VH_VL=(file_repo+"\\all_match_file\\C_file\\VH_VL_file\\all_combined")
    result_folder=(file_repo+"\\result")
    removed_record_folder=(file_repo+"\\all_match_file\\C_file\\C_CDR3")
    #首先处理C开头的CDR3的序列
    for folder_n in [C_CDR3_folder,VH_CDR3_folder,VH_VL,C_VH_VL]:
        if os.path.exists(folder_n):
            print folder_n
            os.chdir(folder_n)
            files_all=os.listdir(folder_n)
            for file_n in files_all:
                if file_n.find("_sequences")!=-1 and file_n.find("freq")!=-1:
                    folder_name=file_n.split("_sequences")[0]
                    if folder_n.find("C_CDR3")!=-1:
                        target_name=(folder_name+"_CDR3_removed_species_freq.fasta")
                        target_path_name=(result_folder+"\\"+folder_name+"\\"+"CDR3_frequency"+"\\"+target_name)
                        
                    elif folder_n.find("VH_CDR3")!=-1:
                        target_name=(folder_name+"_CDR3_freq.fasta")
                        target_path_name=(result_folder+"\\"+folder_name+"\\"+"CDR3_frequency"+"\\"+target_name)

                        
                    elif folder_n.find("combined_VH_VL")!=-1:
                        target_name=(folder_name+"_VH_VL_freq.fasta")
                        target_path_name=(result_folder+"\\"+folder_name+"\\"+"CDR_frequency"+"\\"+target_name)

                        
                    elif folder_n.find("VH_VL_file")!=-1:
                        target_name=(folder_name+"_VH_VL_removed_species_freq.fasta")
                        target_path_name=(result_folder+"\\"+folder_name+"\\"+"CDR_frequency"+"\\"+target_name)
                    if os.path.exists(target_path_name):
                        os.remove(target_path_name)
                    shutil.copyfile(file_n,target_path_name)
    if oa.path.exists(removed_record_folder):
        os.chdir(removed_record_folder)
        shutil.copy("sequence_number_removed.txt",result_folder)
                


#记录给定文件的共有序列，并且返回序列的列表
def get_common_seq(file_path,files_all):#files_all是文件名称的列表
    os.chdir(file_path)
    common_name=""
    set_all_seq=set()
    for file_n in files_all:
        all_seq=[]
        file_name=os.path.splitext(file_n)[0]
        common_name=common_name+"_"+file_name
        with open(file_n,"r") as file_content:
            seq_all=file_content.readlines()
        for seq_n in seq_all:
            if not seq_n.startswith(">"):
                all_seq.append(seq_n)
                #print all_seq
        if not set_all_seq:
            set_all_seq=set(all_seq)
            print set_all_seq
            print len(set_all_seq)
        elif set_all_seq:
            set_all_seq=set_all_seq&set(all_seq)
    #print set_all_seq
    list_common_seq=list(set_all_seq)
    common_seq_file=open("%s_common_seq.txt"%common_name,"w")
    i=0
    for seq_n in list_common_seq:
        i+=1
        common_seq_file.write(str(i)+"\t"+seq_n)
    common_seq_file.close()
    print len(list_common_seq)
    return list_common_seq
            
#定义一个功能找出与common_seq不同的序列并且写到新的文件中
def find_unique_seq(file_path,new_file,common_seq):
    os.chdir(file_path)
    file_name,file_type=os.path.splitext(new_file)
    assert file_type==".fasta","the file should be a fasta format"
    seq_unique=open("%s_unique_sequences.fasta"%file_name,"w")
    with open(new_file,"r") as file_content:
        seq_all=file_content.readlines()
    for i in xrange(len(seq_all)):
        seq_n=seq_all[i]
        if not seq_n.startswith(">"):
            seq_id=seq_all[i-1]
            if common_seq.count(seq_n)==0:
                seq_unique.write(seq_id)
                seq_unique.write(seq_n)
    seq_unique.close()
 
                
def rm_empty_file(file_path):
    os.chdir(file_path)
    files_all=os.listdir(file_path)
    for file_n in files_all:
        if not os.path.getsize(file_n):
            os.remove(file_n)

        
    
    


        
        



#首先第一步解析hmm输出文件，并且缓存序列文件
#首先转入到目标文件夹，并且读取相应的hmmer输出的文件
#下面是自定义文件名和路径名，file_path,hmmer_output,fasta_file分别是hmmer输出文件所在的文件夹，hmmer输出的文件名，含有氨基酸序列的fasta文件名
if __name__=="__main__":
    
    #2018.12.5获取唯一序列
    #file_path=r"E:\Program\project\yanli\FGS_2018.12.5\unique_seq_analysis"#这一行不注释
    #files_all=["Smix1234.fasta","20181017-CD-S-mix2-1-2.fasta"]
    #common_seq=get_common_seq(file_path,files_all)
    #print common_seq
    #files_all=["20181023-P1-CD-ma-ScFv.fasta","20181030-P1-CD-ma-ScFv.fasta","20181106-50nM-CD-ma-Fab.fasta"]
    #common_seq=get_common_seq(file_path,files_all)
    #print common_seq
    #new_file="20181023-P1-CD-ma-ScFv.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #new_file="20181030-P1-CD-ma-ScFv.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #new_file="20181106-50nM-CD-ma-Fab.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #files_all=["20181023-P1-CD-ma-ScFv.fasta","20181030-P1-CD-ma-ScFv.fasta"]
    #common_seq=get_common_seq(file_path,files_all)
    #print common_seq
    #new_file="20181023-P1-CD-ma-ScFv.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #new_file="20181030-P1-CD-ma-ScFv.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    
    #files_all=["20181023-P1-CD-ma-ScFv.fasta","20181106-50nM-CD-ma-Fab.fasta"]
    #common_seq=get_common_seq(file_path,files_all)
    #print common_seq
    #new_file="20181023-P1-CD-ma-ScFv.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #new_file="20181106-50nM-CD-ma-Fab.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #files_all=["20181030-P1-CD-ma-ScFv.fasta","20181106-50nM-CD-ma-Fab.fasta"]
    #common_seq=get_common_seq(file_path,files_all)
    #print common_seq
    #new_file="20181030-P1-CD-ma-ScFv.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #new_file="20181106-50nM-CD-ma-Fab.fasta"
    #find_unique_seq(file_path,new_file,common_seq)
    #file_path=r"E:\Program\project\yanli\FGS_2018.12.5\unique_seq_analysis"#这一行不注释
    #files_all=["test1.fasta","test2.fasta"]
    #common_seq=get_common_seq(file_path,files_all)
    #folder_all=[r"E:\Program\project\yanli\FGS_2018.12.5\unique_seq_analysis\comparison_of_3",r"E:\Program\project\yanli\FGS_2018.12.5\unique_seq_analysis\comparision_of_2\20181023-P1-CD-ma-ScFv_and_20181030-P1-CD-ma-ScFv",r"E:\Program\project\yanli\FGS_2018.12.5\unique_seq_analysis\comparision_of_2\20181023-P1-CD-ma-ScFv_and_20181106-50nM-CD-ma-Fab",r"E:\Program\project\yanli\FGS_2018.12.5\unique_seq_analysis\comparision_of_2\20181030-P1-CD-ma-ScFv_and_20181106-50nM-CD-ma-Fab"]
    #for folder_n in folder_all:
    #    len_out=3
    #    freq_analysis_fasta_ID(folder_n,len_out)
        
    


    
    ##########第九次运行2018.12.5
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
    file_repo=r"E:\Program\project\yanli\FGS_2018.12.6\all_files"#这个文件夹是所有含有测序数据的文件夹所在的文件夹
    #定义每一个含有一代测序文件的文件夹，并且对其进行分析
    #首先定义一个总的目录，在这里所有含有一代测序文件的文件夹都包含在其中
    result_folder=(file_repo+"\\result")
    all_fasta=(result_folder+"\\fasta_file")
    all_folders=get_folder(file_repo)
    print all_folders
    #organize_result(file_repo,all_folders)#该步骤运行完之后及时注释，防止覆盖原来的文件
    scfv_analysis(file_repo,all_folders)

    #########################第二步###########################################
    #接下来进行HMMSCAN分析，并且将所有的比对结果和fasta文件都拷贝到一个文件夹中
    ###########################################################################



    #############################第三步###########################################
    ##############################################################################
    ######接下来是不用修改的部分（将所有需要新建的文件夹在这里建立，之后假如有需要修改的话相应的软件中也需要修改）########################
    #match_repo=(file_repo+"\\all_match_file")#这一行不用注释
    #VH_folder=(match_repo+"\\VH")#这一行不用注释
    #VL_folder=(match_repo+"\\VL")#这一行不用注释
    #hmm_out_folder=(file_repo+"\\all_hmm_out")#这个文件夹需要从其他地方拷贝过来然后解压在file_repo文件夹中，该行不用注释
    #VH_VL_folder=(match_repo+"\\VH_VL_file")#这一行不用注释
    #VH_CDR3=(VH_folder+"\\VH_CDR3")#这个文件夹将会进行CDR3的频率分析
    #C_file=(match_repo+"\\C_file")#这个文件夹将会储存所有CDR以及CDR3的分析
    #C_combined=(C_file+"\\VH_VL_file")#这个文件夹将会用来分析C文件的合并频率
    #combined_folder=(VH_VL_folder+"\\combined_VH_VL")#这一行不用注释
    #rm_folder=(C_combined+"\\rm_folder")#这一行不用注释
    #C_combined_all=(C_combined+"\\all_combined")
    #combined_all=(combined_folder+"\\all_combined")
    #C_CDR3=(C_file+"\\C_CDR3")
    #rm_folder_CDR3=(C_CDR3+"\\without_species")
    #species_CDR3_only_folder=(C_CDR3+"\\CDR3_only")
    #CDR3_only_folder=(VH_CDR3+"\\CDR3_only")
    #species_CDR3_formal_id_folder=(C_CDR3+"\\formal_id")
    #CDR3_formal_id_folder=(VH_CDR3+"\\formal_id")
   # 
    
    #mkdir_folder(match_repo)
    #mkdir_folder(VH_folder)
    #mkdir_folder(VL_folder)
    #mkdir_folder(VH_VL_folder)
    #mkdir_folder(VH_CDR3)
    #mkdir_folder(C_file)
    #mkdir_folder(C_combined)
    #mkdir_folder(combined_folder)
    #mkdir_folder(rm_folder)
    #mkdir_folder(C_combined_all)
    #mkdir_folder(combined_all)
    #mkdir_folder(C_CDR3)
    #mkdir_folder(rm_folder_CDR3)
    #mkdir_folder(species_CDR3_only_folder)
    #mkdir_folder(CDR3_only_folder)#这里运行一遍以后及时注释，否则会多次对文件进行提取CDR3的操作
    #mkdir_folder(species_CDR3_formal_id_folder)
    #mkdir_folder(CDR3_formal_id_folder)
    
    
"""    
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

"""    
    






























    
    
    ######第8次运行2018.11.28 (for yanli)
    #本次将会更新程序，程序将能够在CDR序列的记录文件中包含物种信息，轻重链的信息以及evalue信息，并且增加自动创建文件夹进行自动分析的过程
    #第一步利用FGS模块将核苷酸序列合并。并且拆分为VH VL 两个文件，同时记录下来含有终止密码子的序列
    #提前定义好特异性的参数
    
    #F=["CAGTTTTGGCACAGGCGGCC","CAGTTTTGGCACGGGCGGCC","AAAGAGAGGCCGAGGCGGCC","TTTTGGCACAAGCGGCCGCT","CAGTTTTGGCACACGCGGCC","TAGTTTTGGCACAGGCGGCC","AACGCGAGGCCGAGGCGGCC","CAGTTTTGGCCCAGGCGGCC"]
    #R=["GGGCCGTCGGTGGTCGACCT","GGGCCGTCGGTGGTCGACTG","GGGCCGTCGGTGGCGACTGC","GGGCCGTCGGTGGTCAACTG","AGGCCGTCGGTGGTCGACTG","GGCCGTCGGTGGTCGACCTG","GGCCCGTCGGTGGTCGACCT","GCCTCCACAAGGGCCATCGG","CAGTTTTGGCCCAGGCGGCC","GCCTCCACCAAGGGGCCGTC","GGGCCGTCGGTGGTCGACCT"]
    #split_all=["RAGGGDSGGGGSGGLGGAAA","RAGGGGSGGGGSGGLGGAAA","APGGGGSGGGGSGGLGGAA","APGGGGSGGGGSGRTGGAAA","RAGGGGSGGG","RAGGGSGGGGSGGLGGAAA","RAGGGCSGGGGSGGLGGAAA","RAGGEGSGGGGSGGLGGAAA","RAVGGGSGGGGSGGLGGAAA","RAGRGGSGGGGSGELGGAAA","RAGGGGSGGSGSGGLGGAAA","RAGGGSSGGGGSGGLGGAAA","RAGGGGSEGGGSGGLGGAAA"]#这里需要在循环过程中不断地添加新的linker序列
    #正式进行分析
    #首先新建文件夹来储存所有的轻重链的fasta文件
    #import shutil
    #首先定义将用来储存fasta文件的文佳佳所在位置以及异常序列文件夹所在的位置，这一部分运行完后及时注释掉，以免影响后续的数据分析
    #file_repo=r"E:\Program\project\yanli\FGS_2018.11.27"#将来fasta文件夹所在的文件夹和异常序列文件夹所在的文件夹,这一行不用注释
    #fasta_home=(file_repo+"\\fasta_file")# 这一行不用注释
    #abnormal_home=(file_repo+"\\abnormal_seq")#这一行不用注释
    #if os.path.exists(fasta_home):#判断两个文件夹是否存在，假如不存在那么就建立新的文件夹，假如存在那么就是删除后重新建立
    #    print "ooooook"
    #    shutil.rmtree(fasta_home)
    #if os.path.exists(abnormal_home):#判断两个文件夹是否存在，假如不存在那么就建立新的文件夹，假如存在那么就是删除后重新建立
    #    shutil.rmtree(abnormal_home)
    #os.mkdir(fasta_home)
    #os.mkdir(abnormal_home)
    

    #文件夹处理完毕后再进行文件的处理
    ######第一个文件夹，本文件夹中的类型为fab库的重链测序
    #file_path=r"E:\Program\project\yanli\FGS_2018.11.27\ZQ\C2018110231\20181106-50nM-CD-ma-Fab"
    #folder_name=os.path.split(file_path)[-1]#取出当前所在文件夹的文件夹名称，便于后续组织化的储存数据
    #esul=(abnormal_home+"\\"+folder_name)#这里定义储存当前异常序列的文件夹名称
    #os.mkdir()#创建储存当前文件夹下的异常序列的文件夹
    #F=["TTTTGGCACAAGCGGCCGCT","TTTTGGCCCAAGCGGCCGCT"]#VH的上游引物
    #R=["GGGCCGTCGGTGTTCCCCCT","GCCTCCACCAAGGGGCCGTC","GCCTCCACCAAGGGCCCATC","GCCTCCACCAAGGGGCCATC"]#VH的下游引物
    #chain_type="VH"
    #FGS_analysis(file_path,F,R,chain_type)
    #record_abnormal(file_path,F,R)
    #os.chdir(file_path)
    #shutil.copy("sequence_with_stop_codon.fasta",abnormal_n)
    #shutil.copy("%s_sequences_VH.fasta"%folder_name,fasta_home)
    #shutil.copy("%s_sequences_VL.fasta"%folder_name,fasta_home)
    #第二个文件夹
    #file_path=r"E:\Program\project\yanli\FGS_2018.11.27\ZQ\C2018110231\20181017-CD-S-mix2-1-2"
    #folder_name=os.path.split(file_path)[-1] #取出当前所在文件夹的文件夹名称，便于后续组织化的储存数据
    #abnormal_n=(abnormal_home+"\\"+folder_name)#这里定义储存当前异常序列的文件夹名称
    #os.mkdir(abnormal_n)#创建储存当前文件夹下的异常序列的文件夹
    #rev_c=None
    #chain=None
    #FGS_analysis(file_path,F,R,chain,rev_c,split_all)
    #record_abnormal(file_path,F,R)
    #os.chdir(file_path)
    #shutil.copy("sequence_with_stop_codon.fasta",abnormal_n)
    #if os.path.exists("%s_sequences_VH.fasta"%folder_name):
    #    shutil.copy("%s_sequences_VH.fasta"%folder_name,fasta_home)
    #if os.path.exists("%s_sequences_VL.fasta"%folder_name):
    #    shutil.copy("%s_sequences_VL.fasta"%folder_name,fasta_home)
    #第三个文件夹
    #file_path=r"E:\Program\project\yanli\FGS_2018.11.27\ZQ\C2018110231\20181030-P1-CD-ma-ScFv"
    #folder_name=os.path.split(file_path)[-1] #取出当前所在文件夹的文件夹名称，便于后续组织化的储存数据
    #abnormal_n=(abnormal_home+"\\"+folder_name)#这里定义储存当前异常序列的文件夹名称
    #if os.path.exists(abnormal_n):
    #    shutil.rmtree(abnormal_n)
    #    os.mkdir(abnormal_n)#创建储存当前文件夹下的异常序列的文件夹
    #rev_c=None
    #chain=None
    #FGS_analysis(file_path,F,R,chain,rev_c,split_all)
    #record_abnormal(file_path,F,R)
    #os.chdir(file_path)
    #shutil.copy("sequence_with_stop_codon.fasta",abnormal_n)
    #if os.path.exists("%s_sequences_VH.fasta"%folder_name):
    #    shutil.copy("%s_sequences_VH.fasta"%folder_name,fasta_home)
    #if os.path.exists("%s_sequences_VL.fasta"%folder_name):
    #    shutil.copy("%s_sequences_VL.fasta"%folder_name,fasta_home)
    
    #第二步利用得到的VH和VL序列在linux系统中利用hmmscan进行分析


    #第三步利用match_fasta_txt匹配路径中的fasta和txt的匹配，并且利用匹配到的结果进行分析，得到CDR区域的序列
    #首先将所有的匹配文件聚集在一起放在总的文件夹中，命名为all_match_file          ####################这里每次运行完要记得注释否则拷贝过来的文件会丢失需要重新拷贝
    #match_repo=(file_repo+"\\all_match_file")#这一行不用注释
    #if os.path.exists(match_repo):
    #    shutil.rmtree(match_repo)
    #os.mkdir(match_repo)#创建储存当前文件夹下的异常序列的文件夹
    #随后将所有的hmm输出文件以及fasta文件拷贝到match_repo文件夹中
    #os.chdir(fasta_home)
    #files_all=os.listdir(fasta_home)
    #for file_n in files_all:
    #    shutil.copy(file_n,match_repo)
    #分别新建VH和VL文件夹将重链和轻链的序列分别进行处理
    #首先新建文件夹
    #os.chdir(match_repo)
    #VH_folder=(match_repo+"\\VH")
    #VL_folder=(match_repo+"\\VL")
    #C_folder=(match_repo+"\\C_file")
    #if os.path.exists(VH_folder):
    #    shutil.rmtree(VH_folder)
    #os.mkdir(VH_folder)
    #if os.path.exists(VL_folder):
    #    shutil.rmtree(VL_folder)
    #os.mkdir(VL_folder)
    #if os.path.exists(C_folder):
    #    shutil.rmtree(C_folder)
    #os.mkdir(C_folder)
    #然后将VH和VL文件分别拷贝到VH和VL文件中，分类的标准是是否有VH或者VL在ID中
    #并且将带有C的文件全部考到Cfolder中进行独特的分析
    #files_all=os.listdir(match_repo)
    #for file_n in files_all:
    #    file_type=os.path.splitext(file_n)[1]
    #    if file_type==".fasta" or file_type==".txt":
    #        if file_n.find("VH")!=-1:
    #            shutil.copy(file_n,VH_folder)
    #        elif file_n.find("VL")!=-1:
    #            shutil.copy(file_n,VL_folder)
    
    ######################分别对VH和VL进行分析，首先分析VH
    #VH_folder=(match_repo+"\\VH")#这一行不用注释
    #VL_folder=(match_repo+"\\VL")#这一行不用注释
    #path=VH_folder
    #scheme="Longest_CDR_H"
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)
    #######################然后分析VL
    #path=VL_folder
    #scheme="Longest_CDR_L"
    #dict_match=match_fasta_txt(path)
    #print dict_match
    #get_cdr_from_hmmer(path,dict_match,scheme)
    #为了下一步的进行需要新建文件夹并且将所有的VH和VL的CDR合并
    #VH_VL_folder=(match_repo+"\\VH_VL_file")#这一行不用注释
    #if os.path.exists(VH_VL_folder):
    #    shutil.rmtree(VH_VL_folder)
    #os.mkdir(VH_VL_folder)
    #os.chdir(VH_folder)
    #files_all=os.listdir(VH_folder)
    #for file_n in files_all:
    #    if file_n.find("CDR_sequences")!=-1:
    #        shutil.copy(file_n,VH_VL_folder)
#
    #os.chdir(VL_folder)
    #files_all=os.listdir(VL_folder)
    #for file_n in files_all:
    #    if file_n.find("CDR_sequences")!=-1:
    #        shutil.copy(file_n,VH_VL_folder)
    #单独取出重链的CDR3并且进行多样性的分析，策略为首先将所有的重链的CDR序列取出来放在VH_CDR3文件夹中，随后将C开头的序列进行分解，将其中的不含兔抗的序列提取出来在进行CDR3的分析
    #VH_CDR3=(VH_folder+"\\VH_CDR3")#这个文件夹将会进行CDR3的频率分析
    #if os.path.exists(VH_CDR3):
    #    shutil.rmtree(VH_CDR3)
    #os.mkdir(VH_CDR3)
    #os.chdir(VH_folder)
    #files_all=os.listdir(VH_folder)
    #for file_n in files_all:
    #    if file_n.find("CDR_sequences")!=-1:
    #        shutil.copy(file_n,VH_CDR3)
    #然后对重链的CDR序列文件进行分解，原理是新建文件夹在其中保存所有以C开头的CDR序列，随后在对其进行分解并且记录删除序列的个数
    #C_CDR3=(VH_folder+"\\C_CDR3")#这个文件夹将会进行CDR3的频率分析
    #if os.path.exists(C_CDR3):
    #    shutil.rmtree(C_CDR3)
    #os.mkdir(C_CDR3)
    #os.chdir(VH_folder)
    #files_all=os.listdir(VH_folder)
    #for file_n in files_all:
    #    if file_n.find("CDR_sequences")!=-1 and file_n.startswith("C"):
    #        shutil.copy(file_n,C_CDR3)
    #os.chdir(C_CDR3)
    #接下来对C开头序列进行分析
    #species="rabbit"
    #divide_species_cdr(C_CDR3,species)
    #将C开头文件处理后的不含rabbit的的文件收集到一个文件夹中
    #processed_folder=(C_CDR3+"\\processed_file")
    #if os.path.exists(processed_folder):
    #    shutil.rmtree(processed_folder)
    #os.mkdir(processed_folder)
    #os.chdir(C_CDR3)
    #files_all=os.listdir(C_CDR3)
    #for file_n in files_all:
    #    if file_n.find("without")!=-1:
    #        shutil.copy(file_n,processed_folder)
    #接下来将合并后的CDR文件进行处理
    

    
    #第四步利用merge_VH_VL将VH和VL的序列的CDR区域合并,在这个过程中由于Fab库需要提前进行处理ID，因此先处理ID，随后在进行分析
    #combined_folder=(VH_VL_folder+"\\combined_VH_VL")#这一行不用注释
    #if os.path.exists(combined_folder):
    #    shutil.rmtree(combined_folder)
    #os.mkdir(combined_folder)
    #os.chdir(VH_VL_folder)
    #files_all=os.listdir(VH_VL_folder)
    #print files_all
    #for file_n in files_all:
    #    if file_n.endswith(".fasta"):
    #        shutil.copy(file_n,combined_folder)
    #all_path=combined_folder#这一行不用注释
    #rm_file="T"
    #split_char="."
    #pos_start=0
    #pos_end=1
    #split_id(combined_folder,split_char,pos_start,pos_end,rm_file)
    #merge_VH_VL(all_path)
    #对不含有C的文件进行处理
    #no_c_folder=(combined_folder+"\\no_c")
    #rm_file="T"
    #split_char="_"
    #pos_start=0
    #pos_end=-3
    #split_id(no_c_folder,split_char,pos_start,pos_end,rm_file)
    #merge_VH_VL(no_c_folder)

    
    

    #提取含有兔抗的序列id
    #species="rabbit"
    #dict_file_id=get_species_id(VH_VL_folder,species)
    #print dict_file_id
    #rm_species_id(combined_folder,species,dict_file_id)
    
    #将所有合并后的文件转移到新的文件夹中用来提取和分析不含兔抗的序列
    #分析特殊文件
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.27\cdr3_grep\special"
    #merge_VH_VL(all_path)
    #第五步，利用freq_analysis_fasta_ID统计所有的fasta文件的频率（这一次VH直接用之前分析的去除头部冗余序列的fasta文件，所以不用重新分析）
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.27\cdr3_grep\combined_VH_VL"
    #freq_analysis_fasta_ID(all_path)
    

    
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.15\hmm_out\longest_CDR\combined_VH_VL\edit_id"
    #split_char="T7"
    #pos_start=-1
    #pos_end=None
    #rm_file="T"
    #rm_repeat_seq(all_path,rm_file)#移除重复序列这样便于做进化树
    #split_id(all_path,split_char,pos_start,pos_end,rm_file)
    #phylo_draw(all_path)
    #第六步，根据已有的已经含有各个CDR的fasta文件，进行CDR3提取和频率分析,CDR3_analysis中含有各个文件的三个CDR序列的文件
    #all_path=combined_folder
    #rm_file="F"
    #get_cdr3(all_path,rm_file)
    #cdr3_only_folder=(combined_folder+"\\CDR3_only")
    #
    #if os.path.exists(cdr3_only_folder):
    #    shutil.rmtree(cdr3_only_folder)
    #os.mkdir(cdr3_only_folder)
    #os.chdir(combined_folder)
    #files_all=os.listdir(combined_folder)
    #for file_n in files_all:
    #    if file_n.find("CDR3_only")!=-1:
    #        if file_n.endswith(".fasta"):
    #            shutil.copy(file_n,cdr3_only_folder)
    #len_out=1
    #all_path=cdr3_only_folder
    #freq_analysis_fasta_ID(all_path,len_out)
    #all_path=combined_folder
    #len_out=6
    #freq_analysis_fasta_ID(all_path,len_out)


            


    
    #all_path=no_c_folder
    #rm_file="F"
    #get_cdr3(all_path,rm_file)
    #cdr3_only_folder=(no_c_folder+"\\CDR3_only")
   # 
    #if os.path.exists(cdr3_only_folder):
    #    shutil.rmtree(cdr3_only_folder)
    #os.mkdir(cdr3_only_folder)
    #os.chdir(no_c_folder)
    #files_all=os.listdir(no_c_folder)
    #for file_n in files_all:
    #    if file_n.find("CDR3_only")!=-1:
    #        if file_n.endswith(".fasta"):
    #            shutil.copy(file_n,cdr3_only_folder)
    #len_out=1
    #freq_analysis_fasta_ID(all_path,len_out)
    #all_path=no_c_folder
    #len_out=6
    #freq_analysis_fasta_ID(all_path,len_out)
    
    #第七步，合并所有的重链序列，然后取出CDR3,在进行频率分析
    #all_path=r"F:\all_result\FGS_2018.11.15\result\freq_analysis\CDR3_frequency"
    #combined_fasta(all_path)#合并目标文件夹中所有的fasta文件
    #rm_file="T"
    #get_cdr3(all_path,rm_file)
    #freq_analysis_num_id(all_path)
    #len_out=1
    #freq_analysis_fasta_ID(all_path,len_out)

    
    #第八步，对文件夹中所有的氨基酸序列进行比对分析
    #clustal_all_folder(all_path)
    #第九部分（2018.11.23），对所有序列的CDR3提取出来并且进行频率分析
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.23\longest_cdr\combined_VH_VL\CDR3_analysis"
    #rm_file="T"
    #get_cdr3(all_path,rm_file)
    #freq_analysis_num_id(all_path)
    #freq_analysis_fasta_ID(all_path)

    #合并后的文件进行分析
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.21\first\longest_cdr\combined_VH_VL\CDR3_analysis\combined_S24-1_S24-2_L"
    #rm_file="T"
    #combined_fasta(all_path)
    #get_cdr3(all_path,rm_file)
    #freq_analysis_fasta_ID(all_path)
    
    #all_path=r"E:\Program\project\yanli\FGS_2018.11.21\first\longest_cdr\combined_VH_VL\CDR3_analysis\combined_S24-1_S24-2_H"
    #rm_file="T"
    #combined_fasta(all_path)
    #get_cdr3(all_path,rm_file)
    #freq_analysis_fasta_ID(all_path)
    









     
