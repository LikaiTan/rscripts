library(dplyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(stringr)
library(openxlsx)

# all dates

dt1 <- seq(as.Date("20220420",format="%Y%m%d"), as.Date("20220810",format="%Y%m%d"), by=7) %>% format(format ="%Y%m%d" )
dt2 <- seq(as.Date("20220817",format="%Y%m%d"), as.Date("20230125",format="%Y%m%d"), by=7)%>% format(format ="%Y%m%d" )
dt2 <- c(dt2, '20230129')

dt <- c(dt1, dt2)

# read the data

allcsv_raw <- map(dt, ~ read.csv(
  paste0('https://www.coronavirus.gov.hk/files/death_analysis/Case_fatality_rate/Case_fatality_rate_', .x, '.csv'))
) %>% setNames(dt)

# remove table titles
allcsv <- map(allcsv_raw, ~  
               .x[-(1:4), ] %>%  mutate(Version.date. = str_replace(Version.date., 'Any', 'Other')) )



# add date 
# remove does 1 data
# combine all dataframes into one data frame


dose234 <- map2(allcsv,names(allcsv),  ~ .x[c(3, 9, 10, 13:15, 17:19), c(1, 3:11)] %>%  mutate(todate = .y)%>% 
                  `colnames<-`(c('Vaccine',  'deathrate_3to11' ,   
                                 paste0('deathrate_',seq(10,70,10), 'to', seq(19,79,10) ),
                                 'deathrate_over80', 'ToDate')) %>%  filter(grepl('Dose|Unvacc', Vaccine) )  ) %>% 
  reduce(.f =rbind) 


# 
# write.xlsx(dose234, '/home/big/tanlikai/Hongkongvaccinedata_dose2_3_4.xlsx',asTable = T)


dose234$deathrate_over80 <- as.numeric(dose234$deathrate_over80)*100

ggplot(dose234 , aes(x = ToDate, y = as.numeric(deathrate_over80)*100, color =Vaccine , group = Vaccine))+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90))+ylab('deathrate over80(%)')

ggplot(dose234 , aes(x = ToDate, y = as.numeric(deathrate_70to79)*100, color =Vaccine , group = Vaccine))+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90))

allcsv


dose23_70 <- map2(allcsv,names(allcsv),  ~ .x[c(3, 9, 10, 13:15, 17:19), c(1, 10)] %>%  mutate(todate = .y)  ) %>% reduce(.f =rbind) %>% 
  `colnames<-`(c('Vaccine', 'deathrate_70_to_79', 'ToDate'))


dose3data
dose23_70$deathrate_70_to_79 <- as.numeric(dose23_70$deathrate_70_to_79)*100

ggplot(dose23_70 %>%  filter(Vaccine !=  'Unvaccinated#') ,
       aes(x = ToDate, y = deathrate_70_to_79, color =Vaccine, group = Vaccine ))+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90))+ylab('deathrate_70_to_79(%)')



map2(allcsv,names(allcsv),  ~ .x[c(3, 9, 10, 12:15), c(1, 7)] %>%  mutate(todate = .y)  ) %>% reduce(.f =rbind) %>% 
  `colnames<-`(c('Vaccine', 'deathrate_40_to_49', 'ToDate')) %>%  mutate(deathrate_40_to_49 = as.numeric(deathrate_40_to_49)*100)  %>% 
  ggplot(
         aes(x = ToDate, y = deathrate_40_to_49, color =Vaccine, group = Vaccine ))+
  geom_line()+
  theme(axis.text.x = element_text(angle = 90))+ylab('deathrate_40_to_49(%)')


