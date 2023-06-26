from os import listdir
from os.path import isfile, join
import numpy as np
import re
import pandas as pd
import numpy as np

def data_processing(my_path, taxon_level, presence = True, leaf = False, f5488 = False):
    file_list = [f for f in listdir(my_path) if isfile(join(my_path, f))]
    # Python program to convert .tsv file to .csv file
    # importing re library

    line_list = []
    # reading given tsv file
    final_df = pd.DataFrame()
    for tsv in file_list:
        line_list= []
        with open(my_path + tsv, 'r') as myfile: 
            for line in myfile:
                if taxon_level in line and '@' not in line:
                # Replace every tab with comma
                    fileContent = re.sub("\t", ",", line)
                    fileContent = re.sub("\n", "", fileContent)

                    # Writing into list
                    line_list.append(fileContent.split(','))

        
        df = pd.DataFrame(line_list)
        #print(df)

        df = df.rename(columns={3:'taxon', 4:tsv[-6:-4]})
        df1 = df[['taxon', tsv[-6:-4]]]
        # Create the index
        index_ = df1.taxon
        df1 = df1.drop(['taxon'], axis = 1)

        # Set the index
        df1.index = index_
        df1.head()
        df1 = df1.T
        #print('Final DF:')
        #print(final_df)
        #print('DF1:')
        #print(df1)
        #df1 = df1.replace(np.nan, 0)
        #df1 = df1.replace('NaN', 0)
        if presence == True:
            df1 = df1.astype(float)
            df1[df1 > 0] = 1

        final_df = pd.concat([final_df,df1], axis=0, sort = True)
        final_df = final_df.replace(np.nan,0)
        #print(final_df)
    new_list = []
    for i in file_list:
        if leaf != True and f5488 != True:
            new_list.append(i[4:7])
        elif leaf == True:
            new_list.append(i[4:8])
        else:
            if i[9].isalpha():
                new_list.append(i[7:13])
            else:
                new_list.append(i[7:12])
            
    
    R1_R2_list = ['R1', 'R2']*int(len(new_list)/2)
    final_df.insert(0, 'samples', new_list[:], True)
    final_df.insert(1, 'paired_end', R1_R2_list[:], True)
    final_df.reset_index(inplace = True)
    final_df = final_df.drop(['index'], axis = 1)

    return(final_df)
    # output

