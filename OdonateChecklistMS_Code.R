library(tidyr)
library(countrycode)
library(stringr)
library(plyr)
library(lubridate)
library(data.table)
library(purrr)
library(tidyverse)
library(janitor)
library(dplyr)

####FILES NEEDED#####
#Checklist data (literature and occurrences in separate files)
#Master taxonomy (for harmonizing by valid name for the supplementary files)
#Country neighbors border proportions



###Step 1: Bring in species-country-source combinations data table to convert to ISO3 codes###
####Input literature checklist data includes columns for canonical, country, and source ###
Checklist <- OdonataGlobalChecklist_Supplementary1_081621 #insert lit checklist CSV here
###Input occurrence checklist data is cleaned using Stefan's standardized workflow - contact him for that code####
Occurrence <-FinalOdonata_Occurrences070821 # insert occurrence checklist CSV here
#Remove columns as needed
Checklist <- subset(OdonataGlobalChecklist_Supplementary1_081621, select= -c(Source))
####Step 2: Convert the raw country name (from occurrences or checklists) to ISO3 codes #####
updated <- countrycode(Checklist$RawCountry, "country.name", "iso3c")
Checklist$Country <- c(updated)
#Remove the duplicate canonical-country combinations 
Occurrence <- unique(Occurrence)

##For all remaining steps you only need the columns for 'canonical' and 'country' for each data type

######Step 3: Country Interpolation Script############
##Can either use a file that has combined/bound all literature and occurrence countries (most inclusive way of adding additional species-country combination) or just literature or just occurrence checklist data
#Bind all data sources 
PresenceData <- rbind(Checklist,Occurrence)
#Or just use the literature checklist 
PresenceData <- Checklist
##Add file of the proportions of each country's shared borders
CountryNeighborsISO <- CountryNeighborsISO_allborders
################Baselist of possible missing countries##########################
################Baselist of possible missing countries##########################
####Merge the combinations that it's known in, and then you get the neighbors for that
Base_Comb <- merge(Checklist, CountryNeighborsISO_allborders, by='Country', all.x = TRUE)
###Swap the countries it's known in as the neighbor, then the focal country
Base_Comb <- Base_Comb %>% rename(Country = Neighbor, Neighbor = Country)
##Get rid of extra columns
Base_Clean <- subset(Base_Comb, select= -c(Neighbor,land_border_proportion,total_border_proportion,Source))
##List of the species and all the potential focal countries that it could be in, based on the present neighbors
Base_Clean <- Base_Clean[,c("Canonical","Country")]
##Reorders into alphabetical order of the species
Base_cond <- Base_Clean %>% group_by(Canonical) %>% summarise(Country)
##Unique combinations that it could be in (all known + unknown)
Base_cond <- unique(Base_cond)
##Create list of the known species-country combinations, despite source
Ode_Clean <- subset(Checklist, select = -c(Source))
##Create list of the unique countries each species is known in 
Ode_Clean <- unique(Ode_Clean)
###List of only the countries that each species could be known from (possible countries that each is not already in)
Base_Var <- dplyr::setdiff(Base_cond, Ode_Clean) %>% filter(!is.na(Country))

##########Make a df with the possible combinations of countries to neighbors which are missing###############
###Possible combinations of countries to neighbors 
Country_Var <- subset(Base_Comb, select= -c(land_border_proportion,total_border_proportion))
###Get rid of anything that does not have a neighbor (islands)
Country_Var_Clean <-Country_Var %>% filter(!is.na(Country))
##Reorder as the other list
Country_Var_Clean <- Country_Var_Clean[,c("Canonical","Country", "Neighbor")]
###Group by the name and country, then the neighbors = the countries that it is already known from
Country_Var_Cond <- Country_Var_Clean %>% group_by(Canonical,Country) %>% summarise(Neighbor)
####All the one to one permutations that it could be in 
Country_Var_Cond <- unique(Country_Var_Cond)



################################################################################

###Generate list of valid missing countries based on sum of combination of known present neighbor countries####
####For each country we are interested in, neighbor country border proportion of the focal
CountryNeighbor <- subset(CountryNeighborsISO_allborders, select = -land_border_proportion)

#OPTIONAL##
#####Can also subset to be >50% to limit the default additions######
CountryNeighbor <- CountryNeighbor %>% filter(RoundedProportion<=0.49) ###Set Starting parameters
#######

###Return all the rows from dataframe 1 that have matching value in value 2, if multiple matches, has all possible matches (possible combinations to the values for each of the possible combinations)
Combine1 <- inner_join(Country_Var_Cond, CountryNeighbor, by=c("Country","Neighbor")) 
###Return the potential combinations relative to the unique countries that you are interested in 
Combine <- inner_join(Combine1, Base_Var, by=c("Canonical","Country"))
####Neighbor list is all of the present countries, country is all of the focal countries it could be in, rounded proportion is the final sum of rounded proportions coming from the neighbors
condense <- Combine %>% 
  group_by(Canonical, Country) %>% 
  summarise(Neighbor = paste(unique(Neighbor), collapse = ','), total_border_proportion = sum(total_border_proportion)) %>% filter(total_border_proportion>=0.50)
write.csv(condense,"PossibleNeighbors_param0721_checklist.csv")


########SUPPLEMENTARY FILES#############
####Step 4: Get the unique country-name combinations of checklist data to identify union and delta of data ########
#Get the id numbers, valid name, canonical names from the master taxonomy
setDT(mastertaxonomy)
taxocolumns <- mastertaxonomy %>% select(id,accid,family,validName,Canonical)

#Step 5: Harmonize the literature or occurrence data to the master taxonomy (do both source types separately)
setDT(Occurrence)
Occurrence <- Occurrence %>% select(Canonical,gbifID,Latitude,Longitude,Country,coordinateUncertaintyInMeters,elevation,year,eventdate,DataSource)
Occurrence<- left_join(taxocolumns,Occurrence,by = "Canonical")

#Step 6: Get the relevant spatial columns from the checklist data (separately for literature and occurrence data)
Occurrence <- Occurrence %>% select (validName,Country)
Checklist <- Checklist %>% select (validName,Country)

#Step 7: Get the unique country-name combinations of checklist data
# find and eliminate duplicate combinations
Occurrence <- Occurrence[!duplicated(Occurrence[, c("Country", "validName")]),] 
Checklist <- Checklist[!duplicated(Checklist[, c("Country", "validName")]),] 

#Step 8: Put the checklist countries into strings by valid species
Valid_Checklist <- Checklist[,.('LitCountry'=paste(`Country`,collapse= ",")), by = 'validName']
Valid_Occurrence <- Occurrence[,.('OccurCountry'=paste(`Country`,collapse= ",")), by = 'validName']

#Step 9: Merge these to the full valid species taxonomy to create the supplementary table for valid species
Combinations <- merge(Valid_Checklist,Valid_Occurrence, by = "validName", copy=FALSE,all=TRUE, suffix= c(".x", ".y"))
acceptedtaxo <- MOLMasterTaxonomy[ which(MOLMasterTaxonomy$taxonomic_status=='ACCEPTED'), ]
Supplementary <- merge(acceptedtaxo,Combinations, by = "validName", copy=FALSE,all=TRUE, suffix= c(".x", ".y"))

#Step 10: Split the string of countries per species in order to do union and delta (do for each category you are creating)
setDT(Supplementary)
Supplementary <- Supplementary %>%
  mutate(x = strsplit(LitCountry,","),
         y = strsplit(OccurCountry,",")) %>% 
  mutate(UnionLitOccur = map2(x,y,vecsets::vunion)) %>%
  mutate(LitminusOccurrence = map2(x,y,vecsets::vsetdiff)) %>% 
  mutate(OccurrenceminusLit = map2(y,x,vecsets::vsetdiff))  

#Revector
Supplementary$x <- vapply(Supplementary$x, paste, collapse = ", ", character(1L))
Supplementary$y <- vapply(Supplementary$y, paste, collapse = ", ", character(1L))
Supplementary$UnionLitOccur <- vapply(Supplementary$UnionLitOccur, paste, collapse = ", ", character(1L))
Supplementary$LitminusOccurrence <- vapply(Supplementary$LitminusOccurrence, paste, collapse = ", ", character(1L))
Supplementary$OccurrenceminusLit <- vapply(Supplementary$OccurrenceminusLit, paste, collapse = ", ", character(1L))

#Rename all the fields for the supplementary table
names(Supplementary)[names(Supplementary) == "x"] <- "ChecklistCountry"
names(Supplementary)[names(Supplementary) == "y"] <- "OccurrenceCountry"
names(Supplementary)[names(Supplementary) == "UnionLitOccur"] <- "Lit+Occur"
names(Supplementary)[names(Supplementary) == "LitminusOccurrence"] <- "Lit-Occur"
names(Supplementary)[names(Supplementary) == "OccurrenceminusLit"] <- "Occur-Lit"


####Step 11: For each species, make maps of present countries by data source
##Bind the checklist data together and add a column with the source information (source = literature, occurrence, or interpolation)
##You should have columns for species, country, and source (this is to assign the colors of the map)

serrata <- InterpolationSource0721_allpresences %>% filter(validName == "Aeshna serrata Hagen, 1856")

serrata_map <- joinCountryData2Map(serrata, 
                                   joinCode = "ISO3",
                                   nameJoinColumn = "Country")
par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i")
serrata_map1 <- mapCountryData(serrata_map, nameColumnToPlot="Source", 
                               catMethod = "categorical", 
                               missingCountryCol = gray(.8),
                               
                               colourPalette = c("#fc8d59", "#ffffbf", "#91bfdb"), 
                               mapTitle = "Aeshna serrata Hagen, 1856", borderCol = "black")


####Calculate overlap between data sources
#Create dataframes by source
interpolation<-OdonataGlobalChecklist_Supplementary_Supplementary5[(OdonataGlobalChecklist_Supplementary_Supplementary5$Source=="interpolation"),]
literature<-OdonataGlobalChecklist_Supplementary_Supplementary5[(OdonataGlobalChecklist_Supplementary_Supplementary5$Source=="literature"),]
occurrence<-OdonataGlobalChecklist_Supplementary_Supplementary5[(OdonataGlobalChecklist_Supplementary_Supplementary5$Source=="occurrence"),]

#Remove the column of the source
interpolation <-interpolation[, c('validName', 'Country')]
literature <-literature[, c('validName', 'Country')]
occurrence <-occurrence[, c('validName', 'Country')]

#Get only unique rows for the combined literature + occurrence data
lit_and_occur <-rbind(literature,occurrence)
lit_and_occur <- unique(lit_and_occur)

#Species-country in common between literature and occurrence
common_lit_occur <- generics::intersect(literature, occurrence)

#Species-country in literature and not occurrence
lit_not_occur <- anti_join(literature,occurrence)

#Species-country in occurrence and not literature
occur_not_lit <- anti_join(occurrence,literature)

#Number of species in literature
lit_species <-literature[, c('validName')]
lit_species <-unique(lit_species)

#Number of species in occurrences
occur_species <-occurrence[, c('validName')]
occur_species<-unique(occur_species)

#Speciesin common between literature and occurrence
common_sp_lit_occur <- generics::intersect(lit_species, occur_species)

#Get only unique rows for the combined literature + occurrence data
lit_and_occur_sp <-rbind(lit_species,occur_species)
lit_and_occur_sp <- unique(lit_and_occur_sp)




