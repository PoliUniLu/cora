import string 

COLUMN_LABELS = list(string.ascii_uppercase) + ["AA", "BB", "CC", "DD", "EE", "FF"]

def concatenate_strings(arr):
    return ','.join([str(x) for x in arr])

def inclusion_score(arr):
    return sum(arr)/len(arr)

def preprocess_data(df, out_cols,input_cols=None, n_cut=1, inc_score1=1, 
                    case_col=None, inc_score2=None, U=None, inverse_output=False,rename_columns=False):
  
    df_tmp=df.copy()
    if input_cols is None:
        input_cols=[x for x in list(df.columns) if (x not in out_cols) and (x!=case_col)]
        print(input_cols)
    
    
    if(case_col is None or case_col == '-None-'):
        df_tmp['case_col']=df_tmp.index.values
        case_col = 'case_col'
    if (inverse_output):
        df_tmp[out_cols]=abs(df_tmp[out_cols]-1)
    #input_columns = [x for x in df_new.columns if x not in OutCols and x != CaseColl]  
    params = {'inc_score_{}'.format(c) : (c,inclusion_score) for c in out_cols}
    df_grouped=df_tmp.groupby(input_cols).agg(n_cut=(case_col, 'count'), 
                                                 case_ids=(case_col, concatenate_strings),
                                                 **params)
    df_grouped=df_grouped[df_grouped['n_cut']>=n_cut]
    inc_columns=['inc_score_{}'.format(i) for i in out_cols]
    if(inc_score2 is None):
        #print('cicik')
        df_grouped[out_cols]=(df_grouped[inc_columns]>=inc_score1).astype(int)
    else:
        if(U is None):
            raise Exception('When inc.score2 is specified, U must be specified as well.')
        if(U != 0 and U != 1):   
            raise Exception('U must be 0 or 1.')
        if (U == 1):
            df_grouped[out_cols]=(df_grouped[inc_columns]> inc_score2).astype(int)
        if (U == 0):   
            df_grouped[out_cols]=(df_grouped[inc_columns]>= inc_score1).astype(int)
            
        
    res = df_grouped.reset_index()
    
    if rename_columns:
        rename_dic = {k:v for k,v in zip(input_cols, COLUMN_LABELS[:len(input_cols)])}
        res.columns = map(lambda x: rename_dic[x] if x in rename_dic.keys() else x, res.columns)
        print('renamed to : {}'.format(res.columns) )
        print(rename_dic)
        
    return res 
