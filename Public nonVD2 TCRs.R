
# Project: nonVD2 TCR sharing.  -------------------------------------------


dir.SR2 <- '/home/big/tanlikai/Alina_data/Likai_Filespublic/'
SRall<- list.files( path = dir.SR2, pattern = '.txt')[-1]
SRall
#get donor ID 
namesSR <- 
  str_extract(SRall,
              '[:upper:][:upper:][:digit:]{1,4}|Pax[:digit:][:digit:]|primal[:digit:]{1,2}|PK-[:digit:][:digit:]|Pk6') 
namesSR <- namesSR %>% toupper() %>% str_remove('-')
namesSR %>% unique()
# [1] "LT1142" "LT1165" "LT1268" "LT1294" "LT304"  "LT469"  "LT581"  "LT711" 
# [9] "LT781"  "LT801"  "LT915"  "LT937"  "PAX01"  "PAX03"  "PAX04"  "PAX05" 
# [17] "PAX06"  "PAX07"  "PAX08"  "PAX09"  "PAX10"  "PAX11"  "PAX14"  "PAX17" 
# [25] "PAX18"  "PAX19"  "PAX20"  "PAX21"  "PAX22"  "PK14"   "PK21"   "PK24"  
# [33] "PK17"   "PK23"   "PK26"   "PK31"   "PK32"   "PK16"   "PK7"    "PK29"  
# [41] "PK28"   "PK20"   "PK30"   "PK5"    "PK3"    "PK4"    "PK19"   "PK1"   
# [49] "PK22"   "PK2"    "PK8"    "PK15"   "PK10"   "PK11"   "PK27"


###group unique cdr3 based on donor ID
SR_detail_all <- map(SRall, function(x){
  read.table(file = paste0(dir.SR2, '/',x), header = T)  %>%  mutate(txt = x)
}) %>% set_names(SRall) %>% 
  purrr::reduce(.f = rbind) %>% 
  mutate(donor = str_extract(txt, 'Pk6|[:upper:][:upper:][:digit:]{1,4}|Pax[:digit:][:digit:]|primal[:digit:]{1,2}|PK-[:digit:][:digit:]') ) %>%
  mutate(donor = toupper(donor)) %>% mutate(donor = str_remove(donor, '-')) %>% ungroup() %>% 
  
  group_by(cdr3aa, donor,v,d,j) %>% summarise(counts = sum(count)   )  %>%
  group_by(donor) %>% group_split() 


SR_detail_all


SRdonors <- map(SR_detail_all, function(x) unique(x$donor)) %>% as.list()
SR_detail_all <- SR_detail_all %>% set_names(SRdonors) 

SR_detail_all


SR_detail_all_cdr3 <- map(SR_detail_all, ~ 
                            .x %>% mutate(cdr3  =  case_when(
                              v == 'TRDV2' ~ str_replace(cdr3aa, '^\\w\\w',""),
                              !is.na(cdr3aa) ~ as.vector(cdr3aa)
                            ))%>%  pull(cdr3) %>%  unique()
)


length(SRdonors)


SR_detail_all_onlyTCR <- SR_detail_all %>% purrr::reduce(.f = rbind) %>%  pull(cdr3aa) %>% unique()

SR_detail_all_onlyTCR


DVdata_D  <- list.files('/home/big/tanlikai/Fetal_thymus_TCR_David_V/') %>% 
  str_subset('_D.txt')%>% 
  map( ~ read.table(paste0('/home/big/tanlikai/Fetal_thymus_TCR_David_V/', .x), header = T) %>%  
         mutate(sample = .x, 
                donor = str_extract(.x, '[FP]T[123]')) %>%   group_by(cdr3aa, donor,v,d,j) %>% 
         summarise(counts = sum(count)   ) ) 



SRTCRlist <- SR_detail_all %>% map(~ .x$cdr3aa %>% unique) %>%   unlist() 

DVTCRlist <- DVdata_D  %>% map(~ .x$cdr3aa %>% unique) %>%   unlist() 

TCRlist_total <- c(SRTCRlist,DVTCRlist)%>% table() %>% sort(decreasing = T) %>%  as.data.frame() %>% `colnames<-`(c('cdr3aa', 'Public_count'))


DVdata_D  %<>%  map(~ .x %>% left_join(TCRlist_total, by = 'cdr3aa', suffix = c('','') ) %>%  
                      mutate(publicity = case_when(Public_count == 1 ~ 'private',
                                                   nr(Public_count, 2, 5) ~ 'low public (2 ~ 5)', 
                                                   nr(Public_count, 6, 15)  ~ 'moderate public (6 ~ 15)',
                                                   Public_count > 15 ~ 'high public (> 15)')) )


DVdata_D %>%  map(~ .x %>% ungroup() %>% count(donor, v,j, publicity) %>% filter(v== 'TRDV1') %>% arrange(publicity) )


DVdata_D %>% purrr::reduce(.f = rbind) %>%  pull(cdr3aa) %>% table() %>% sort(decreasing = T)

DVdata_D %>% map(~ .x %>%filter(cdr3aa == 'CACDTGSSWDTRQMFF'))

DVdata_D$cdr3aa %>% unique()







# based on nt -------------------------------------------------------------





SR_detail_allnt <- map(SRall, function(x){
  read.table(file = paste0(dir.SR2, '/',x), header = T)  %>%  mutate(txt = x)
}) %>% set_names(SRall) %>% 
  purrr::reduce(.f = rbind) %>% 
  mutate(donor = str_extract(txt, 'Pk6|[:upper:][:upper:][:digit:]{1,4}|Pax[:digit:][:digit:]|primal[:digit:]{1,2}|PK-[:digit:][:digit:]') ) %>%
  mutate(donor = toupper(donor)) %>% mutate(donor = str_remove(donor, '-')) %>% ungroup() %>% 
  
  group_by(cdr3nt, donor,v,d,j) %>% summarise(counts = sum(count)   )  %>%
  group_by(donor) %>% group_split() 

# SR_detail_allnt <- SR_detail_all %>% purrr::reduce(.f = rbind) %>%  pull(cdr3nt) %>% unique()


SR_detail_allnt %>% purrr::reduce(.f = rbind) %>%  pull(cdr3nt) %>% table() %>% sort(decreasing = T)


SR_detail_all_onlyTCRnt <- SR_detail_allnt %>% purrr::reduce(.f = rbind) %>%  pull(cdr3nt) %>% unique()

DVdata_Dnt  <- list.files('/home/big/tanlikai/Fetal_thymus_TCR_David_V/') %>% 
  str_subset('_D.txt')%>% 
  map( ~ read.table(paste0('/home/big/tanlikai/Fetal_thymus_TCR_David_V/', .x), header = T) %>%  
         mutate(sample = .x, 
                donor = str_extract(.x, '[FP]T[123]')) %>%   group_by(cdr3nt, donor,v,d,j) %>% 
         summarise(counts = sum(count)   ) ) 

DVdata_Dnt %>% purrr::reduce(.f = rbind) %>%  pull(cdr3nt) %>% table() %>% sort(decreasing = T) %>% max()



publicnt <- DVdata_Dnt %>%map(~ .x %>%  filter(cdr3nt %in% SR_detail_all_onlyTCRnt) %>% ungroup %>%  count(v) )
publicnt
