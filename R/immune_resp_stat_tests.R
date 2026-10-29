# Statistical comparisons of mucosal immunity (IgA binding and pseudo-neutralisation):
#   1. Children vs infants at each state / sample / time / poliovirus type
#   2. Baseline vs follow-up within each age group
# Results are written to output/compare_age_groups.xlsx and output/compare_time.xlsx

source(here::here("R", "libraries.R"))

load(here(path("data","anon_neut.Rdata")))
load(here(path("data","anon_vacc.Rdata")))
load(here(path("data","anon_bind.Rdata")))

# Add age group (Infants/Children) to the neutralisation data
df_neut<-df_neut|>
  left_join(df_vacc|>select(age_cat,
                            id),
            by="id")
# ---- 1. Compare age groups ----

# Pseudo-neutralisation: number and % positive in each age group, and Fisher's
# exact test comparing the two groups. Done for the primary cutoff (no suffix)
# and the sensitivity-analysis cutoff (suffix "_SA").
neut<-
  df_neut|>
  pivot_longer(c(neut_pos,neut_pos_SA),names_to="Cutoff")|>
  mutate(Cutoff=ifelse(Cutoff=="neut_pos","","_SA"))|>
  group_by(state_print,age_cat,sample_print,time_print,pv_print,Cutoff)|>
  summarise(pos=sum(value,na.rm = T),
            n=sum(!is.na(value)))|>
  pivot_wider(names_from=age_cat,values_from=c(pos,n))|>
  group_by(pv_print,sample_print,state_print,time_print,Cutoff)|>
  mutate(
    p_compare_age=fisher.test(matrix(c(pos_Children,n_Children-pos_Children,pos_Infants,n_Infants-pos_Infants),ncol=2))$p.value,
    Children=paste0(pos_Children,"/",n_Children," (",signif(pos_Children/n_Children*100,2),"%)"),
    Infants=paste0(pos_Infants,"/",n_Infants," (",signif(pos_Infants/n_Infants*100,2),"%)"))|>
  select(-starts_with("pos"),-starts_with("n"))|>
  pivot_wider(values_from=c(p_compare_age,Children,Infants),names_from=Cutoff,names_sep="")

# Add age group to the IgA binding data
df_bind<-
  df_bind|>
  left_join(df_vacc|>select(age_cat,
                            id),
            by="id")

# Welch two-sample t-test of IgA binding (bind_trans), children vs infants,
# within one state / sample type / time point / poliovirus type. Returns the p-value.
age_compare_t<-
  function(data,st,sa,ti,pv){
    t.test(data$bind_trans[data$age_cat=="Children"&
                             data$state_print==st&
                             data$sample_print==sa&
                             data$time_print==ti&
                             data$pv_print==pv],
           data$bind_trans[data$age_cat=="Infants"&
                             data$state_print==st&
                             data$sample_print==sa&
                             data$time_print==ti&
                             data$pv_print==pv],paired=F)$p.value
  }

# One row per combination to test, with the t-test p-value
frm<-
  distinct(df_bind,state_print,sample_print,time_print,pv_print)
frm$p_compare_age<-sapply(1:nrow(frm),function(i){age_compare_t(data=df_bind,
                                                                frm$state_print[i],
                                                                frm$sample_print[i],
                                                                frm$time_print[i],
                                                                frm$pv_print[i])})

# Add median (IQR) and n of IgA binding for each age group
frm<-
  frm|>
  left_join(
    df_bind|>
      group_by(state_print,
               sample_print,
               time_print,
               pv_print,
               age_cat)|>
      filter(!is.na(bind_trans))|>
      summarise(med=median(bind_trans),
                p25=quantile(bind_trans,0.25),
                p75=quantile(bind_trans,0.75),
                n=n(),
                data=paste0(round(med,2)," (",round(p25,2),", ",round(p75,2),") n=",n))|>
      select(-med,-p25,-p75,-n)|>
      pivot_wider(names_from=age_cat,values_from=data))

# Stack the IgA and pseudo-neutralisation results and save
cbind(frm,type="IgA",test="t-test",population="All")|>
  bind_rows(cbind(neut,type="Pseudo-neutralization",test="Fishers exact test",population="All"))|>
  arrange(pv_print)|>
  write_xlsx(here(path("output","compare_age_groups.xlsx")))



# ---- 2. Compare time points (baseline vs follow-up) ----

# Pseudo-neutralisation: McNemar's test on paired positive/negative status
# (children with both baseline "b" and follow-up "f" results only), for the
# primary cutoff and the sensitivity-analysis cutoff. Also adds the n/N (%)
# positive at each time point.
compare_neut_time<-function(df_neut){
  bind_rows(
    
    df_neut|>
      select(-time_print,-neut,-neut_pos_SA,-log2_neut)|>
      mutate(neut_pos=factor(neut_pos,levels=c(TRUE,FALSE),labels=c("Pos","Neg")))|>
      pivot_wider(names_from=time,values_from=neut_pos)|>
      filter(!is.na(b),
             !is.na(f))|>
      group_by(pv_print,age_cat,sample_print,state_print)|>
      summarise(p_compare_time=mcnemar.test(b,f)$p.value,
                n=n(),
                test="McNemar's test"),
    
    df_neut|>
      select(-time_print,-neut,-neut_pos,-log2_neut)|>
      mutate(neut_pos_SA=factor(neut_pos_SA,levels=c(TRUE,FALSE),labels=c("Pos","Neg")))|>
      pivot_wider(names_from=time,values_from=neut_pos_SA)|>
      filter(!is.na(b),
             !is.na(f))|>
      group_by(pv_print,age_cat,sample_print,state_print)|>
      summarise(p_compare_time=mcnemar.test(b,f)$p.value,
                n=n(),
                test="McNemar's test (Cut-off sensitivity analysis)"),
  ) |>
    left_join(
      df_neut|>
        pivot_longer(c(neut_pos,neut_pos_SA),names_to="Cutoff")|>
        mutate(test=ifelse(Cutoff=="neut_pos","McNemar's test","McNemar's test (Cut-off sensitivity analysis)"))|>
        group_by(state_print,age_cat,sample_print,time_print,pv_print,test)|>
        summarise(pos=sum(value,na.rm = T),
                  n=sum(!is.na(value)),
                  print=paste0(pos,"/",n," (",signif(pos/n*100,2),"%)"))|>
        select(-pos,-n)|>
        pivot_wider(names_from=time_print,values_from=print)
    )
}


# IgA binding: paired t-test of baseline vs follow-up (children with both
# results only), plus median (IQR) and n at each time point.
compare_bind_time<-function(df_bind){
  df_bind|>
    select(id,time,bind_trans,pv_print,age_cat,sample_print,state_print)|>
    pivot_wider(names_from=time,values_from=bind_trans)|>
    filter(!is.na(b),
           !is.na(f))|>
    group_by(pv_print,age_cat,sample_print,state_print)|>
    summarise(p_compare_time=t.test(b,f,paired=T)$p.value,
              n=n(),
              test="Paired t-test")|>
    left_join(
      df_bind|>
        group_by(state_print,
                 sample_print,
                 time_print,
                 pv_print,
                 age_cat)|>
        filter(!is.na(bind_trans))|>
        summarise(med=median(bind_trans),
                  p25=quantile(bind_trans,0.25),
                  p75=quantile(bind_trans,0.75),
                  n=n(),
                  data=paste0(round(med,2)," (",round(p25,2),", ",round(p75,2),") n=",n))|>
        select(-med,-p25,-p75,-n)|>
        pivot_wider(names_from=time_print,values_from=data))}

# Stack the IgA and pseudo-neutralisation results and save
time<-cbind(compare_bind_time(df_bind),type="IgA",population="All")|>
  bind_rows(cbind(compare_neut_time(df_neut),type="Pseudo-neutralization",population="All"))|>
  arrange(pv_print)

time|>
  write_xlsx(here(path("output","compare_time.xlsx")))
