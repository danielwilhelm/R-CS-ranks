library(dplyr)
math_scores <- readxl::read_excel(file.path("data-raw", "pisa2022.xls"), 
                                  na=c("", "—", "†"), 
                                  sheet = "Report 1-Table Cleared") %>% 
  rename(math_score="Average", math_se="Standard Error", Year="Year/Study")
reading_scores <- readxl::read_excel(file.path("data-raw", "pisa2022.xls"), 
                                     na=c("", "—", "†"), 
                                     sheet = "Report 2-Table Cleared") %>%
  rename(reading_score="Average", reading_se="Standard Error",
         Year="Year/Study")
science_scores <- readxl::read_excel(file.path("data-raw", "pisa2022.xls"),  
                                     na=c("", "—", "†"), 
                                     sheet = "Report 3-Table Cleared") %>%
  rename(science_score="Average", science_se="Standard Error",
         Year="Year/Study")

pisa_2022_2018 <- dplyr::full_join(science_scores, reading_scores) %>%
  dplyr::full_join(math_scores)

pisa2022 <- dplyr::filter(pisa_2022_2018, Year == 2022) %>%
  dplyr::select(-Year)
usethis::use_data(pisa2022)