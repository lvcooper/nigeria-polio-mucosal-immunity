# Figure 1: (a) map of study states and LGAs; (b) timeline of nOPV2 campaigns
# (SIAs), cohort birth dates and sample collection, with the prevalence of type 2
# stool neutralising activity at each survey

source(here::here("R", "libraries.R"))

# ---- Event timeline data ----------------------------------------------------

# One row per event (SIA round or survey) per state, with start and end dates.
# Label each row as an SIA or a survey, and name SIAs by the vaccine used.
df<-read_xlsx(here(path("data","events.xlsx")))|>
  mutate(start=as.Date(start),
         end=as.Date(end),
         event_type=ifelse(grepl("SIA",event),"SIA","Survey"),
         event=ifelse(grepl("SIA",event),vaccine,event))

# Right-hand end of the cohort bars
lim<-as.Date("2024-03-01")

# Add the birth-date range of each cohort, derived from age at baseline.
# For each cohort the "Oldest" and "Youngest" rows hold the earliest and latest
# possible birth dates (baseline date minus age in months):
#   Children 12-24 months at baseline, Infants 6-10 months at baseline.
# These rows are collapsed into a single range per cohort in df_birth below.
df<-
  bind_rows(df,
            df|>
              filter(event=="Baseline assessment")|>
              mutate(end=lim,
                     start=start-24*365/12+1,
                     event="Oldest",
                     event_type="Child cohort"),
            df|>
              filter(event=="Baseline assessment")|>
              mutate(end=lim,
                     start=start-12*365/12,
                     event="Youngest",
                     event_type="Child cohort"),
            df|>
              filter(event=="Baseline assessment")|>
              mutate(end=lim,
                     start=start-10*365/12+1,
                     event="Oldest",
                     event_type="Infant cohort"),
            df|>
              filter(event=="Baseline assessment")|>
              mutate(end=lim,
                     start=start-6*365/12,
                     event="Youngest",
                     event_type="Infant cohort"))


load(here(path("data","anon_neut.Rdata")))
load(here(path("data","anon_vacc.Rdata")))

# Proportion of stool samples with type 2 neutralising activity (neut_pos) by
# state, age group and survey, with exact binomial 95% CIs. Each point is placed
# at the midpoint of its survey's dates.
df_points<-
  df_neut|>
  left_join(df_vacc|>
              select(id,age_cat,bopv_ri_3plus,ri_info_source))|>
  filter(pv==2,sample=="s")|>
  summarise(.by=c(pv_print,sample_print,state_print,age_cat,time_print),
            pos=sum(neut_pos,na.rm=T),
            n=sum(!is.na(neut_pos)),
            p=pos/n,
            lci=binom.test(pos,n)$conf[1],
            uci=binom.test(pos,n)$conf[2])|>
  # Attach the survey date (time_print matches "Baseline"/"Follow-up")
  left_join(df|>
              filter(event_type=="Survey")|>
              mutate(.by=c(event,state),
                     time_print=event|>str_remove(" assessment"),
                     state_print=state,
                     date=mean(c(start,end)))|>
              select(time_print,state_print,date))|>
  mutate(event=age_cat)


# SIA rounds, drawn as arrows at the midpoint of each campaign
df_sia<-
  df|>
  filter(event_type=="SIA")|>
  mutate(.by=start,
         state_print=state,
         event=paste(event,event_type),
         date=mean(c(start,end)))

# Cohort birth-date ranges (one bar per cohort and state) plus the survey
# sample-collection windows. event_type_num sets each bar's vertical position.
df_birth<-
  df|>
  filter(grepl("cohort",event_type))|>
  # Collapse the Oldest/Youngest rows into one range per cohort
  summarise(.by=c(state,event_type),
         end=max(start),
         start=min(start))|>
  mutate(.by=start,
         date=mean(c(start,end)))|>
  bind_rows(df|>
              filter(grepl("Survey",event_type)))|>
  mutate(state_print=state,
         event=recode(event_type,
                      "Child cohort"="Dates of birth (Children)",
                      "Infant cohort"="Dates of birth (Infants)",
                      "Survey"="Sample collection"),
         event_type_num=as.factor(event_type)|>as.numeric())

# Panel b: one facet per state. Bars show cohort birth dates and sample
# collection, arrows show SIAs, and points with CIs show the prevalence of
# type 2 neutralising activity at each survey.
# Bars, arrows and labels sit at fixed heights (y about 0.25-0.45) above the data
# points; these y values are positioning offsets, not data.
timeline_alt<-
  ggplot()+
  geom_segment(data=df_sia,aes(x=date,xend=date,y=0.3,yend=0.25,color=event),arrow=arrow(length=unit(0.1,"inches")))+
  geom_segment(data=df_birth,aes(x=start,xend=end,y=0.4-event_type_num*0.02,yend=0.4-event_type_num*0.02,color=event),linewidth=5)+
  # Survey rows have no `date`, so geom_text drops them (the warning about
  # 4 removed rows); the "Sample collection" label is added by annotate() below
  geom_text(data=df_birth,aes(x=date,y=0.45-event_type_num*0.02,label=event),size=3)+
  annotate("text",x=as.Date("2023-12-01"),y=0.45-3*0.02,label="Sample collection",size=3)+
  scale_x_date()+
  scale_alpha_manual(values=c(1,0.5))+
  scale_color_manual(values=brewer.pal(6,"Dark2")[c(2,3,2,3,4,5,6)])+
  theme(legend.position="top")+
  guides(alpha="none")+
  labs(x="",y="Prevalence of type 2 mucosal neutralising activity (titres >=8)",color="",shape="",linetype="")+
  geom_pointrange(data=df_points,aes(x=date,y=p,ymin=lci,ymax=uci,color=event,shape=event,linetype=event),position=position_dodge(width=14))+
  facet_wrap(~state_print,ncol=1)+theme_bw()

# ---- Map data ---------------------------------------------------------------

# Colours: pal[1]/pal[2] for Sokoto (study LGAs / rest of state),
# pal[3]/pal[4] for Adamawa
pal<-c(brewer.pal(3,"Dark2"),brewer.pal(3,"Set2"))[c(1,4,2,5)]

# GADM 4.1 Nigeria boundaries (https://gadm.org; https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_NGA.gpkg)
# ADM_0 = country, ADM_1 = states, ADM_2 = LGAs
gadm_file<-here(path("data","gadm41_NGA.gpkg"))
adm0<-st_read(gadm_file,layer="ADM_ADM_0",quiet=TRUE)
adm1<-st_read(gadm_file,layer="ADM_ADM_1",quiet=TRUE)
# Keep LGAs in the two study states and compute label positions (centroids)
adm2<-st_read(gadm_file,layer="ADM_ADM_2",quiet=TRUE)|>
  filter(NAME_1%in%c("Sokoto","Adamawa"))|>
  mutate(CENTER_LON=st_coordinates(st_centroid(geom))[,1],
         CENTER_LAT=st_coordinates(st_centroid(geom))[,2])

# Study LGAs
sokoto_lgas<-c("Sokoto North","Sokoto South","Wamako")
yola_lgas<-c("Yola North","Yola South")

# Zoom windows for the insets: bounding box of the study LGAs plus 0.2 degrees
box1<-st_bbox(adm2|>filter(NAME_2%in%sokoto_lgas))+0.2*c(-1,-1,1,1)
box2<-st_bbox(adm2|>filter(NAME_2%in%yola_lgas))+0.2*c(-1,-1,1,1)

# Panel a: Nigeria with the two study states shaded and the study LGAs
# highlighted. The invisible point just fixes the plot extent.
main_map<-ggplot()+
  geom_point(data=cbind.data.frame(x=1,y=6),aes(x=x,y=y),alpha=0)+
  geom_sf(data=adm1,fill="white")+
  geom_sf(data=adm0,fill='transparent')+
  geom_sf(data=adm2|>filter(NAME_1=="Sokoto"),fill=pal[2],alpha=0.5)+
  geom_sf(data=adm2|>filter(NAME_1=="Adamawa"),fill=pal[4],alpha=0.5)+
  geom_sf(data=adm2|>filter(NAME_2%in%sokoto_lgas),fill=pal[1])+
  geom_sf(data=adm2|>filter(NAME_2%in%yola_lgas),fill=pal[3])+
  theme_void()

# Insets: the main map zoomed to each state's study LGAs, with LGA names
inset1<-main_map+
  geom_text_repel(data=adm2|>filter(NAME_2%in%sokoto_lgas)|>arrange(NAME_2),
                  aes(x=CENTER_LON,y=CENTER_LAT,label=NAME_2),
                  # Order: Sokoto North, Sokoto South, Wamako. Nudged apart so
                  # the leader lines don't cross
                  nudge_x=c(0.12,0.12,-0.1),nudge_y=c(0.1,-0.1,0.1),
                  min.segment.length=0,force=0.5,size=3)+
  coord_sf(xlim=box1[c(1,3)],ylim=box1[c(2,4)])+
  theme(panel.border=element_rect(fill=NA))+
  labs(subtitle="Sokoto")
inset2<-main_map+
  geom_text_repel(data=adm2|>filter(NAME_2%in%yola_lgas),
                  aes(x=CENTER_LON,y=CENTER_LAT,label=NAME_2),min.segment.length=0.5,force=500,size=3)+
  coord_sf(xlim=box2[c(1,3)],ylim=box2[c(2,4)])+
  theme(panel.border=element_rect(fill=NA))+
  labs(subtitle="Adamawa")

# ---- Combine and save -------------------------------------------------------

# Place the two insets on top of the main map (x, y, width, height are
# fractions of the panel)
w<-0.4
map<-ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset1, x = 0.05, y = 0.6, width = w, height = w)+
  draw_plot(inset2, x = 0.55, y = 0.05, width = w, height = w)+
  theme(plot.background = element_rect(fill="white",color=NA))

# Map (a) above timeline (b); ggsave saves the last plot drawn
ggarrange(map,timeline_alt,ncol=1,heights=c(0.4,0.6),labels="auto")
ggsave(filename=here(path("output","figure1.png")),width=8,height=8,dpi=500)
ggsave(filename=here(path("output","figure1.pdf")),width=8,height=8)
