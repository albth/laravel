import pandas as pd, numpy as np
def load():
    c=pd.read_excel('ctx.xlsx',sheet_name='|FC| >1')
    c=c.rename(columns={'condition-JD130_CtxvsJD1884_Ctx.log2FC':'lfc','condition-JD130_CtxvsJD1884_Ctx.FDR':'fdr','condition-JD130_CtxvsJD1884_Ctx.RawP':'p'})
    c['lfc']=-c['lfc']   # flip: positive = vcAst (JD1884) enriched
    pan=[f'JD130_Ctx_{i}' for i in range(1,5)]; vc=[f'JD1884_Ctx_{i}' for i in range(1,5)]
    c['pan']=c[pan].mean(1); c['vc']=c[vc].mean(1)
    c=c[['Gene','lfc','p','fdr','pan','vc']].copy(); c['region']='Ctx'
    s=pd.read_excel('str.xlsx',sheet_name='|FC| >2')
    s=s.rename(columns={'condition-JD130_StrvsJD1884_Str.log2FC':'lfc','condition-JD130_StrvsJD1884_Str.FDR':'fdr','condition-JD130_StrvsJD1884_Str.RawP':'p'})
    s['lfc']=-s['lfc']
    pan=[f'JD130_Str_{i}.norm' for i in range(1,5)]; vc=[f'JD1884_Str_{i}.norm' for i in range(1,5)]
    s['pan']=s[pan].mean(1); s['vc']=s[vc].mean(1)
    s=s[['Gene','lfc','p','fdr','pan','vc']].copy(); s['region']='Str'
    return c,s
