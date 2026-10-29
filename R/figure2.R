# Figure 2: stool IgA binding vs stool neutralisation titre, by poliovirus type

source(here::here("R", "libraries.R"))

load(here(path("data","anon_neut.Rdata")))
load(here(path("data","anon_vacc.Rdata")))
load(here(path("data","anon_bind.Rdata")))

# One row per child x poliovirus type x sample x time, with IgA binding
# (bind_trans) alongside the neutralisation titre (log2_neut) and age group
df_comp<-
  df_bind|>
  left_join(df_vacc|>select(id,age_cat))|>
  left_join(df_neut)|>
  select(id,state_print,time_print,pv_print,sample_print,bind_trans,log2_neut,age_cat)|>
  mutate(log2_neut=as.numeric(log2_neut))

# Correlation between IgA and neutralisation for each poliovirus type and
# sample type (stool/nasal), pooling states, time points and age groups.
df_cor<-
  df_comp|>
  group_by(pv_print,sample_print)|>
  summarise(N=n(),
            p=cor.test(bind_trans,log2_neut,method='spearman',exact=F)$p.value,
            rho=cor.test(bind_trans,log2_neut)$est,
            p_num=p,
            p=ifelse(p<0.0001,"p<0.0001",paste0("p=",round(p,3))))

# Stool samples with a neutralisation result only; n is the number of samples
# at each titre within each poliovirus type (used for the boxplot filter and
# the "n=" labels)
df<-
  df_comp|>
  filter(sample_print=="Stool")|>
  filter(!is.na(log2_neut))|>
  mutate(.by=c(log2_neut,pv_print),n=n())

p<-ggplot()+
  # Boxplots only where a titre has more than 10 samples
  geom_boxplot(data=df|>
                 filter(n>10),aes(y=bind_trans,x=log2_neut,group=log2_neut),outlier.alpha=0)+
  # All individual samples, jittered horizontally because titres are discrete
  geom_jitter(data=df,aes(y=bind_trans,x=log2_neut),shape=1,height=0,width=0.25,alpha=0.5)+
  # Correlation p-value in each panel; rho is shown only if p<0.05
  geom_text(data=df_cor|>
              filter(sample_print=="Stool")|>
              mutate(est=case_when(p_num<0.05~paste0(", rho=",round(rho,2)),
                                   p_num>=0.05~""),
                     lab=paste0(p,est)),
            aes(x=7,y=0,label=lab),size=3)+
  # Sample count at each titre
  geom_text(data=df|>
    count(pv_print,log2_neut)|>
              mutate(lab=paste0("n=",n)),
            aes(x=log2_neut,y=5,label=lab),size=2.5)+
  # Dashed line between log2 titres 2 and 3 (titres 4 and 8), i.e. the
  # seropositivity cut-off of >=8
  geom_vline(data=df|>distinct(pv_print)|>mutate(xin=2.5),aes(xintercept=xin),alpha=0.5,linetype=2)+
  facet_grid(~pv_print)+
  # x-axis is on the log2 scale but labelled with the actual titre (2^x)
  scale_x_continuous(breaks=1:10,labels=2^(1:10),minor_breaks=1:10)+
  labs(x="Stool neutralisation titre",y="Log10 stool IgA adjusted MFI")+
  theme_bw()

ggsave(filename=here(path("output","figure2.pdf")),p,width=9,height=3)
