# Tables of immune responses (neutralising titres and IgA) by state, age group
# and poliovirus type, written as Word documents to output/

source(here::here("R", "libraries.R"))

load(here(path("data","anon_neut.Rdata")))
load(here(path("data","anon_bind.Rdata")))
load(here(path("data","anon_vacc.Rdata")))

# Wide table with one row per child x poliovirus type x sample type, holding
# baseline (_b) and follow-up (_f) neutralisation and IgA (_iga) results plus
# vaccination history. Columns derived at the end:
#   conv   = seroconversion: negative at baseline and positive at follow-up
#            (NA if already positive at baseline)
#   boost  = >=4-fold rise in neutralising titre from baseline to follow-up
#   either = seroconversion or 4-fold rise
table_base<-
  df_neut|>
  pivot_wider(names_from=time,values_from=c(neut_pos,neut_pos_SA,neut),id_cols=c(pv,id,sample_print))|>
  full_join(df_bind|>
              pivot_wider(names_from=time,values_from=bind_trans,id_cols=c(pv,id,sample_print),names_glue="{time}_iga"),
            by=c("id","sample_print","pv"))|>
  left_join(df_vacc|>select(age_cat,state_print,nopv2_sia,bopv_ri_3plus_num,ipv_1plus_num,ri_info_source,id),by="id")|>
  mutate(cohort=paste0(state_print," ",age_cat),
         conv=case_when(!neut_pos_b&neut_pos_f~T,
                        !neut_pos_b&!neut_pos_f~F,
                        neut_pos_b~NA,
                        is.na(neut_f)|is.na(neut_b)~NA),
         boost=case_when((neut_f/neut_b)>=4~T,
                         (neut_f/neut_b)<4~F,
                         is.na(neut_f)|is.na(neut_b)~NA),
         either=case_when(conv==T~T,
                          boost==T~T,
                          conv==F&is.na(boost)~F,
                          conv==F&boost==F~F,
                          boost==F&is.na(conv)~F,
                          is.na(neut_f)|is.na(neut_b)~NA))

# Table labels: the primary seropositivity cutoff and the sensitivity-analysis
# cutoff (_SA) are read from the data; labels are keyed by column name (the
# _1/_2/_3 suffixes are poliovirus types)
neut_cutoff<-unique(df_neut$cutoff)
neut_cutoff_SA<-unique(df_neut$cutoff_SA)
labs<-list(neut_pos_b=paste0("Neutralising titres >=",neut_cutoff," at baseline"),
           neut_pos_f=paste0("Neutralising titres >=",neut_cutoff," at follow-up"),
           neut_pos_SA_b=paste0("Neutralising titres >=",neut_cutoff_SA," at baseline"),
           neut_pos_SA_f=paste0("Neutralising titres >=",neut_cutoff_SA," at follow-up"),
           conv=paste0("<",neut_cutoff," to >=",neut_cutoff," neutralising titres baseline to follow-up"),
           boost="4-fold rise in neutralising titres baseline to follow-up",
          either=paste0("Either <",neut_cutoff," to >=",neut_cutoff," or 4-fold rise in neutralising titres"),
                                            b_iga="IgA log10 MFI at baseline",
                                            f_iga="IgA log10 MFI at follow-up",
          neut_pos_b_1=paste0("Type 1 neutralising titres >=",neut_cutoff," at baseline"),
          neut_pos_f_1=paste0("Type 1 neutralising titres >=",neut_cutoff," at follow-up"),
          neut_pos_SA_b_1=paste0("Type 1 neutralising titres >=",neut_cutoff_SA," at baseline"),
          neut_pos_SA_f_1=paste0("Type 1 neutralising titres >=",neut_cutoff_SA," at follow-up"),
          b_iga_1="Type 1 IgA log10 MFI at baseline",
          f_iga_1="Type 1 IgA log10 MFI at follow-up",
          neut_pos_b_2=paste0("Type 2 neutralising titres >=",neut_cutoff," at baseline"),
          neut_pos_f_2=paste0("Type 2 neutralising titres >=",neut_cutoff," at follow-up"),
          neut_pos_SA_b_2=paste0("Type 2 neutralising titres >=",neut_cutoff_SA," at baseline"),
          neut_pos_SA_f_2=paste0("Type 2 neutralising titres >=",neut_cutoff_SA," at follow-up"),
          b_iga_2="Type 2 IgA log10 MFI at baseline",
          f_iga_2="Type 2 IgA log10 MFI at follow-up",
          neut_pos_b_3=paste0("Type 3 neutralising titres >=",neut_cutoff," at baseline"),
          neut_pos_f_3=paste0("Type 3 neutralising titres >=",neut_cutoff," at follow-up"),
          neut_pos_SA_b_3=paste0("Type 3 neutralising titres >=",neut_cutoff_SA," at baseline"),
          neut_pos_SA_f_3=paste0("Type 3 neutralising titres >=",neut_cutoff_SA," at follow-up"),
          b_iga_3="Type 3 IgA log10 MFI at baseline",
          f_iga_3="Type 3 IgA log10 MFI at follow-up")

# Table 2: Type 2 response by age and state before and after
# (columns = sample type x state x age group; continuous IgA shown as median (IQR))
table_base|>
  filter(pv==2)|>
  mutate(cohort=paste0(sample_print," ",state_print," ",age_cat))|>
  select(cohort,neut_pos_b,neut_pos_f,neut_pos_SA_b,neut_pos_SA_f,b_iga,f_iga)|>
  tbl_summary(by=cohort,
                           label=labs,
                           statistic = list(
                             all_continuous() ~ "{median} ({p25}, {p75}) n={N_nonmiss}",
                             all_categorical() ~ "{n} / {N} ({p}%)",
                             all_dichotomous() ~ "{n} / {N} ({p}%)"
                           ),
              digits = list("{N_nonmiss}"=0,"{median}"=1,"{p25}"=1,"{p75}"=1))|># remove the missing row for the test variables
  remove_row_type(type = "missing")%>%
  as_gt()%>%
  gtsave(filename=here(path("output","table2_type2_age_state.docx")))

# Type 2 responses stratified by age, state and nopv2_sia status (supplementary
# table 3); nopv2_sia is whether the child reportedly received the most recent
# nOPV2 SIA dose
table_base|>
  filter(pv==2)|>
  mutate(cohort=paste0(sample_print," ",state_print," ",age_cat," ",nopv2_sia))|>
  select(cohort,neut_pos_b,neut_pos_f,neut_pos_SA_b,neut_pos_SA_f,b_iga,f_iga)|>
  tbl_summary(by=cohort,
                           label=labs,
                           statistic = list(
                             all_continuous() ~ "{median} ({p25}, {p75}) n={N_nonmiss}",
                             all_categorical() ~ "{n} / {N} ({p}%)",
                             all_dichotomous() ~ "{n} / {N} ({p}%)"
                           ),
              digits = list("{N_nonmiss}"=0,"{median}"=1,"{p25}"=1,"{p75}"=1))|># remove the missing row for the test variables
  remove_row_type(type = "missing")%>%
  as_gt()%>%
  gtsave(filename=here(path("output","table2_type2_age_state_nopv2.docx")))

# Type 1 & 3 response in fully vaccinated children: 3+ bOPV and 1+ IPV doses
# documented on a vaccination card. Types 1 and 3 are shown side by side by
# widening the data to one row per child.
table_base|>
  filter(pv!=2,
         bopv_ri_3plus_num==1,
         ipv_1plus_num==1,
         ri_info_source=="Card")|>
  mutate(cohort=paste0(sample_print," ",state_print," ",age_cat))|>
  select(id,cohort,pv,neut_pos_b,neut_pos_f,neut_pos_SA_b,neut_pos_SA_f,b_iga,f_iga)|>
  pivot_wider(names_from=pv,values_from=c(neut_pos_b,neut_pos_f,neut_pos_SA_b,neut_pos_SA_f,b_iga,f_iga))|>
  select(cohort,neut_pos_b_1,neut_pos_f_1,neut_pos_SA_b_1,neut_pos_SA_f_1,b_iga_1,f_iga_1,
         neut_pos_b_3,neut_pos_f_3,neut_pos_SA_b_3,neut_pos_SA_f_3,b_iga_3,f_iga_3)|>
  tbl_summary(by=cohort,
              label=labs,
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75}) n={N_nonmiss}",
                all_categorical() ~ "{n} / {N} ({p}%)",
                all_dichotomous() ~ "{n} / {N} ({p}%)"
              ),
              digits = list("{N_nonmiss}"=0,"{median}"=1,"{p25}"=1,"{p75}"=1))|># remove the missing row for the test variables
  remove_row_type(type = "missing")%>%
  as_gt()%>%
  gtsave(filename=here(path("output","type13_age_state_ri_card.docx")))

# Type 1 and 3 response overall population (all children, regardless of
# vaccination history)
table_base|>
  filter(pv!=2)|>
  mutate(cohort=paste0(sample_print," ",state_print," ",age_cat))|>
  select(id,cohort,pv,neut_pos_b,neut_pos_f,neut_pos_SA_b,neut_pos_SA_f,b_iga,f_iga)|>
  pivot_wider(names_from=pv,values_from=c(neut_pos_b,neut_pos_f,neut_pos_SA_b,neut_pos_SA_f,b_iga,f_iga))|>
  select(cohort,neut_pos_b_1,neut_pos_f_1,neut_pos_SA_b_1,neut_pos_SA_f_1,b_iga_1,f_iga_1,
         neut_pos_b_3,neut_pos_f_3,neut_pos_SA_b_3,neut_pos_SA_f_3,b_iga_3,f_iga_3)|>
  tbl_summary(by=cohort,
              label=labs,
              statistic = list(
                all_continuous() ~ "{median} ({p25}, {p75}) [n={N_nonmiss}]",
                all_categorical() ~ "{n} / {N_nonmiss} ({p}%)",
                all_dichotomous() ~ "{n} / {N_nonmiss} ({p}%)"
              ),
              digits = list("{N_nonmiss}"=0,"{median}"=1,"{p25}"=1,"{p75}"=1))|># remove the missing row for the test variables
  remove_row_type(type = "missing")%>%
  as_gt()%>%
  gtsave(filename=here(path("output","type13_age_state_full_pop.docx")))

