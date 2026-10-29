# Participant characteristics and vaccination history tables, by state and age
# group. Output: output/table1.docx and output/tableS4.docx

source(here::here("R", "libraries.R"))

load(here(path("data","anon_vacc.Rdata")))
load(here(path("data","anon_neut.Rdata")))

# Prepare variables for the summary tables
df<-
  df_vacc|>
  # Collapse parental education into two levels
  mutate(mothers_education=mothers_education|>
           factor(levels=c("Never attended",
                           "Non-formal",
                           "Primary",
                           "Secondary",
                           "Tertiary",
                           "Postgraduate"),
                  labels=c("Primary or below",
                           "Primary or below",
                           "Primary or below",
                           "Secondary or above",
                           "Secondary or above",
                           "Secondary or above")),
         fathers_education=fathers_education|>
           factor(levels=c("Never attended",
                           "Non-formal",
                           "Primary",
                           "Secondary",
                           "Tertiary",
                           "Postgraduate"),
                  labels=c("Primary or below",
                           "Primary or below",
                           "Primary or below",
                           "Secondary or above",
                           "Secondary or above",
                           "Secondary or above")),
         # Split routine immunisation history by source (parental recall vs
         # vaccination card) so each is summarised separately
         bopv_recall=case_when(ri_info_source=="Recall"~bopv_ri_3plus),
         bopv_card=case_when(ri_info_source=="Card"~bopv_ri_3plus),
         ipv_recall=case_when(ri_info_source=="Recall"~ipv_1plus),
         ipv_card=case_when(ri_info_source=="Card"~ipv_1plus),
        # Order age groups so Infants appear before Children in the table
        age_cat=factor(age_cat,levels=c("Infants","Children"),ordered=T))|>
  # id is kept so the supplementary table can be restricted to children with
  # neutralisation results; it is dropped before each table is built
  select(id,
sex,
num_child_5_yrs_in_hh,
mothers_education,
fathers_education,
ri_info_source,
bopv_recall,
bopv_card,
ipv_recall,
ipv_card,
sia_doses_FU_cat,
nopv2_sia,
fipv_sia_recall,
state_print,
age_cat)|>
  # Row labels shown in the tables
  set_variable_labels(sex="Sex (female)",
                      num_child_5_yrs_in_hh="Number children under five in household",
                      mothers_education="Mother's education (primary or below)",
                      fathers_education="Father's education (primary or below)",
                      bopv_recall="bOPV history (recall, 3 or more)",
                      bopv_card="bOPV history (card, 3 or more)",
                      ipv_recall="IPV history (recall, 1 or more)",
                      ipv_card="IPV history (card, 1 or more)",
                      nopv2_sia="Most recent nOPV2 SIA dose received at follow-up (recall)",
                      sia_doses_FU_cat="OPV SIA doses received in lifetime at follow-up (recall)",
                      fipv_sia_recall="Most recent fIPV SIA dose received at follow-up (recall)",
                      ri_info_source="Routine immunization information source")

# Table 1: one block per state, columns by age group. Binary variables show
# only the named level (e.g. "Female", "3+" doses) as n / N (%); the missing
# row is dropped for the vaccination history variables.
df|>
  select(-id)|>
  tbl_strata(strata=state_print,
             .tbl_fun=
               ~.x|>tbl_summary(
                 by = age_cat,
                 type = list(num_child_5_yrs_in_hh ~ "continuous",
                             sex ~ "dichotomous",
                             ipv_recall ~ "dichotomous",
                             bopv_recall ~ "dichotomous",
                             ipv_card ~ "dichotomous",
                             bopv_card ~ "dichotomous",
                             mothers_education ~ "dichotomous",
                             fathers_education ~ "dichotomous"),
                 value = list(sex="Female",
                              mothers_education="Primary or below",
                              fathers_education="Primary or below",
                              bopv_card="3+",
                              bopv_recall="3+",
                              ipv_card="1+",
                              ipv_recall="1+"),
                 statistic = list(
                   num_child_5_yrs_in_hh ~ "{median} ({min}, {max})",
                   all_categorical() ~ "{n} / {N} ({p}%)",
                   all_dichotomous() ~ "{n} / {N} ({p}%)"
                 ),
                 digits = list(all_continuous() ~ 0)
               ) |>
               remove_row_type(variables = c(num_child_5_yrs_in_hh,
                                             bopv_card, 
                                             bopv_recall, 
                                             ipv_card, 
                                             ipv_recall), type = "missing")
  )|>
  as_gt()|>
  gtsave(filename=here(path("output","table1.docx")))
 

 # Reported nOPV2 SIA coverage per LGA: n and % of children by response
  df_vacc |> 
    group_by(lga, nopv2_sia) |> 
    summarise(Freq = n(), .groups = "drop") |> 
    group_by(lga) |> 
    mutate(Proportion = round((Freq*100) / sum(Freq), 3))

# Participant characteristics for whom mucosal neutralising activity was measured.
df|>
  filter(id%in%df_neut$id)|>
  select(-id)|>
  tbl_strata(strata=state_print,
             .tbl_fun=
               ~.x|>tbl_summary(
                 by = age_cat,
                 type = list(num_child_5_yrs_in_hh ~ "continuous",
                             sex ~ "dichotomous",
                             ipv_recall ~ "dichotomous",
                             bopv_recall ~ "dichotomous",
                             ipv_card ~ "dichotomous",
                             bopv_card ~ "dichotomous",
                             mothers_education ~ "dichotomous",
                             fathers_education ~ "dichotomous"),
                 value = list(sex="Female",
                              mothers_education="Primary or below",
                              fathers_education="Primary or below",
                              bopv_card="3+",
                              bopv_recall="3+",
                              ipv_card="1+",
                              ipv_recall="1+"),
                 statistic = list(
                   all_continuous() ~ "{median} ({min}, {max})",
                   all_categorical() ~ "{n} / {N} ({p}%)",
                   all_dichotomous() ~ "{n} / {N} ({p}%)"
                 ),
                 digits = all_continuous() ~ 0
               ) |># remove the missing row for the test variables
               remove_row_type(variables = c(num_child_5_yrs_in_hh,
                                             bopv_card, 
                                             bopv_recall, 
                                             ipv_card, 
                                             ipv_recall), type = "missing")
  )%>%
  as_gt()%>%

  gtsave(filename=here(path("output","tableS4.docx")))

# Supplementary table 5 - cross tab of parent reported nOPV2 doses in focal SIA and in lifetime up to follow-up:
df_vacc %>%
select(nopv2_sia, sia_doses_FU_cat) %>%
tbl_summary(
  by = nopv2_sia,
  percent = "column"
) %>%
add_overall()
