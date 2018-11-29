# 2018.11.19 Summarize Bismark output
# author: Ziao Lin
import sys
import pandas as pd
import numpy as np

read_name_list = []
strand_list = []
chrom_list = []
pos_list = []
methy_list = []

in_file = sys.argv[1]
sample_id = sys.argv[2]

#print len(sys.argv)

##with open("sample_bam/CpG_OT_CW3_TP2_RRBS.txt") as infile:
with open(in_file) as infile:

     next(infile)
     for line in infile:
         a = line
         b = a.split()
         read_name_list.append(b[0])
         strand_list.append(b[1])
         chrom_list.append(b[2])
         pos_list.append(b[3])
         methy_list.append(b[4])

read_name_df = pd.DataFrame(read_name_list)
strand_df = pd.DataFrame(strand_list)
chrom_df = pd.DataFrame(chrom_list)
pos_df = pd.DataFrame(pos_list)
methy_df = pd.DataFrame(methy_list)
combine_df = pd.concat([read_name_df, strand_df, chrom_df, pos_df, methy_df], axis =1)
combine_df.columns = ["read_name", "strand", "chrom", "pos", "methylated"]

#it is too slow for using the set.
#unique_name = list(set(combine_df['read_name']))

#it is too slow to subset dataframe every time.

#instead, since the text is ordered by read names, we can use a counter to check the one that's different from the previous one.
new_read_name_list = []
new_strand_list = []
new_chrom_list = []
new_pos_list = []
new_methy_list = []


previous_read_name = 0
combine_pos = str()
combine_methy = str()

for i in range(len(combine_df)):
    #print (i)
    if combine_df["read_name"].iloc[i] == previous_read_name:
        current_pos = combine_df['pos'].iloc[i]
        combine_pos = combine_pos + str(current_pos) + "\t"
        current_methy = combine_df['methylated'].iloc[i]
        combine_methy = combine_methy + str(current_methy) + "\t"
        
    else:
        previous_read_name = combine_df["read_name"].iloc[i]
        new_pos_list.append(combine_pos)
        new_methy_list.append(combine_methy)
        
        combine_pos = str()
        combine_methy = str()
        
        #no subsetting this time.
        #current_df = combine_df.loc[combine_df["read_name"] == i]
        #current_df = combine_df.loc[combine_df["read_name"] == combine_df["read_name"].iloc[i]]
        
        current_read_name = combine_df['read_name'].iloc[i]
        current_strand = combine_df['strand'].iloc[i]
        current_chrom = combine_df['chrom'].iloc[i]
        current_pos = combine_df['pos'].iloc[i]
        combine_pos = combine_pos + str(current_pos) + "\t"
        
        current_methy = combine_df['methylated'].iloc[i]
        combine_methy = combine_methy + str(current_methy) + "\t"
        
        #current_strand = current_df['strand'].iloc[0]
        #current_chrom = current_df['chrom'].iloc[0]
        
        #current_pos = np.array(current_df['pos'])
    
        #combine_pos = str()
        #for j in current_pos:
            #combine_pos = combine_pos + str(j) + "\t"
        
        #current_methy = np.array(current_df['methylated'])
    
        #combine_methy = str()
        #for k in current_methy:
            #combine_methy = combine_methy + str(k) + "\t"
        
        new_read_name_list.append(current_read_name)
        new_strand_list.append(current_strand)
        new_chrom_list.append(current_chrom)
    

#should save the last combine_pos and combine_methy.
new_pos_list.append(combine_pos)
new_methy_list.append(combine_methy)

new_pos_list = new_pos_list[1:]
new_methy_list = new_methy_list[1:]

new_read_name_df = pd.DataFrame(new_read_name_list)
new_strand_df = pd.DataFrame(new_strand_list)
new_chrom_df = pd.DataFrame(new_chrom_list)
new_pos_df = pd.DataFrame(new_pos_list)
new_methy_df = pd.DataFrame(new_methy_list)


new_combine_df = pd.concat([new_read_name_df, new_strand_df, new_chrom_df, new_pos_df, new_methy_df], axis = 1)
new_combine_df.columns = ["read_name", "strand", "chrom", "pos", "methy_status"]

#write output dataframe.
new_combine_df.to_csv(str(sample_id) + str("_pdr.txt"), sep='\t', index = False)



