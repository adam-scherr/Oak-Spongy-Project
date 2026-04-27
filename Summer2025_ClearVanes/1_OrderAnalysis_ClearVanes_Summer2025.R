###1_OrderAnalysis_ClearVanes_Summer2025
###Order level diversity analysis of Summer 2025 clear vane traps
###April 8, 2026 (but based on clear vane traps from summer 2025)
###Adam Scherr

###Make Data Wide####
##load in your dataset
long.data <- read.csv("~/Downloads/6_Arthropod ID - Clear Vane Sample Order IDs.csv")
##cleaning data so only beetle families are used
library("dplyr")
library("tidyverse")

#check for typos to make sure orders are spelled correctly
attach(long.data)
unique(Order) 
detach(long.data)

#remove blank rows and unknowns
long.data.1 <- long.data %>%
  filter(Order != "", Order != "Unknown")

#check if only the orders of interest are present
attach(long.data.1)
unique(Order) 
detach(long.data.1)

  #it appears that you misspelled Staphylinidae (Staphylindae), Pyrochroidae (Pyrchroidae),
  #Laemophloeidae (Laemophleoidae), Cryptophagidae (Cryptohagidae), Scarabaeidae (Scarabeidae)
  #fix this now in your original excel sheet...okay, now it's been fixed

#remove the unknowns
long.data.2 <- long.data.1%>%
  filter(Order != "Unknown")
attach(long.data.2)
unique(Order) #all unknown families are removed (juveniles, partial specimens); 57 families total
unique (Forest)
unique(BurnType)
unique(TimePoint)
detach(long.data.2)

##make the data frame wide; give each order its own column
str(long.data.2) #Count is listed as a character because of the thrips you did not count. you need it as an integer

#set Count as an integer
long.data.2$Count <- as.integer(long.data.2$Count)
str(long.data.2)

#now make the dat frame wider
wide.data <- long.data.2 %>%
  group_by(Forest, BurnType, TakedownDate, TimePoint, TreeNumber, TrapDestroyed, Order)%>%
  summarize(count = sum(Count))%>% #make new variable called "count" based on sums of Count column
  pivot_wider(names_from = Order, #get names from Family variable
              names_sort = TRUE, #have columns form in alphabetical order
              values_fn = sum, #this line of code does nothing, I tested it
              values_from = count, #define the value for the Family columns
              values_fill = 0) %>% #turn blanks into zeros
  ungroup()
    #23 columns total. First 5 columns are labels. Remaining 18 columns are orders. Looks good.


#confirming that the "values_fn = sum" line of code does nothing.
wide.data.test <- long.data.2 %>%
  group_by(Forest, BurnType, TakedownDate, TimePoint, TreeNumber, TrapDestroyed, Order)%>%
  summarize(count = sum(Count))%>%
  pivot_wider(names_from = Order, #get names from Family variable
              names_sort = TRUE, #have columns form in alphabetical order
              values_fn = sum, 
              values_from = count,
              values_fill = 0) %>%
  ungroup()


identical (wide.data, wide.data.test) #yup, that line of code makes no difference here

##our dataset is ready to be used for analysis! Let's just write it into a fresh CSV to upload to github
getwd()
 #working directory is currently my Clear Vane ID_Summer2025 folder. That will do for now.

# library("readr")
# write_csv(wide.data, "wide.order.data.csv")

###Diversity Analyses####
#remove unusable data points, namely if traps were destroyed
wide.data.1 <- wide.data %>%
  filter(TrapDestroyed == 'no')

#for easier referencing in graphs, make combined forest-site-tree#-timepoint
wide.data.1<- wide.data.1 %>%
  mutate(Unique = paste(Forest, BurnType, TimePoint, TreeNumber, sep = '-'))

#preliminary plot of the data
long.data.3 <- long.data.2 %>%
  filter(Order != "Thysanoptera") #Thysanoptera data is unusable anyway

ggplot(long.data.3, aes(x = Order, y = Count))+
  geom_point()+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
  #this graph doesn't add the counts. It gives you an idea of the data.
  #You can see that Acari, Coleoptera, Collembola, Diptera, Hemiptera, and Hymenoptera
    #are the most abundant orders

###Shannon Diversity####
##load in the vegan package with the Shannon diversity index function
library(vegan)

##make shannon diversity dataset
shan_data <- wide.data.1 [, 7:23] #select just the data containing insect order data; excluding Thysanoptera data
shan_index <- diversity(shan_data, index = "shannon") #run shannon index

#add shannon index values to wide.data.1, make a new dataframe
data_shan <- wide.data.1 %>%
  add_column(shan = shan_index)

#change variable classes
data_shan$Forest <- as.factor(data_shan$Forest)
data_shan$BurnType <- as.factor(data_shan$BurnType)
class(data_shan$Forest)
class(data_shan$BurnType)

#change burn type appearance order
data_shan$BurnType <- factor(data_shan$BurnType, 
                             levels = c("recent", "medium", "unburned"))

#make Shannon diversity index model
library(glmmTMB)
library(lme4)
library(AICcmodavg)
shan_model<-glmmTMB(shan ~ BurnType + (1|Forest), data = data_shan, 
                    family = "gaussian") #(1|Forest) means that Forest is a random effect. BurnType is a fixed effect.
#shan_glm <- lmer(shan ~ BurnType + (1 | Forest), data = data_shan)
shan_model.gaus<-glmmTMB(shan ~ BurnType + (1|Forest), data = data_shan, 
                         family = "gaussian")
shan_model.nb1 <- update(shan_model.gaus, family = "nbinom1")
shan_model.nb2 <- update(shan_model.gaus, family = "nbinom2")
shan_model.nb1.zi <- update(shan_model.nb1, ziformula = ~1) 
shan.model.nb2.zi <- update(shan_model.nb2, ziformula=~1)
shan_model.pois <- update(shan_model.gaus, family = "poisson") 
shan_model.pois.zi <-  update(shan_model.gaus, family = "poisson", ziformula=~1)

list <- list(shan_model.gaus = shan_model.gaus, shan_model.nb1=shan_model.nb1,
             shan_model.nb2=shan_model.nb2, shan_model.nb1.zi=shan_model.nb1.zi,
             shan_model.pois=shan_model.pois, shan_model.pois.zi=shan_model.pois.zi,
             shan.model.nb2.zi = shan.model.nb2.zi)
aictab(list) #lowest aic was for shan_model.gaus
shan_glm<- shan_model.gaus

##test all the residuals
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = shan_glm)
plot(simulationOutput)
  #Everything looks good, but we will test each feature individually to be sure

#test the residual features individually:
testUniformity(simulationOutput) # p-value not significant, residuals are normal
testDispersion(simulationOutput) # p-value not significant, residuals are not overdispersed
testZeroInflation(simulationOutput) #p-value not significant, data not zero inflated
testOutliers(simulationOutput) #now outliers marked in red, no outliers found

##run an ANOVA on your Shannon diversity index glm
library(car)
#?Anova()
anova_shan<-Anova(shan_glm, type = "III") #use type III ANOVA if the interaction 
#of your two independent variables (Order and BurnType) is significant.
#I am still not sure if the interaction is significant or not. You may want to use type = 'II'
anova_shan
summary(anova_shan)

#make summary data containing just the means and standard errors for count data
library(tidyverse)
library(plotrix)

shan_summary <- data_shan %>%
  dplyr::group_by(BurnType) %>%
  dplyr::summarise(
    mean = mean(shan, na.rm = TRUE),
    se = std.error(shan, na.rm = TRUE),
    .groups = "keep"
  )

# Get estimated marginal means (emmeans) for the Treatment factor
library(emmeans)
library(multcomp)
library(multcompView)

shan_emmeans <- emmeans(shan_glm, specs = pairwise ~ BurnType,
                        adjust = "sidak")
shan_cld_results <- cld(shan_emmeans, adjust = "sidak", Letters = letters)
shan_letters_df <- as.data.frame(shan_cld_results)
print(names(shan_letters_df))

shan_letters_df <- shan_letters_df %>%
  dplyr::select(BurnType, .group) %>%
  dplyr::rename(Letter = .group)

# Merge letters with summary data
shan_summary <- shan_summary %>%
  left_join(shan_letters_df, by = "BurnType")
View(shan_summary)

shan_summary$BurnType <- factor(shan_summary$BurnType, 
                                levels = c("recent", "medium", "unburned"))
shan_graph <- shan_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 6, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE)+
  ggtitle("Shannon Diversity Index")
shan_graph
  #At the order level, the Shannon Diversity did not differ based on burn treatment

###Count Data Diversity####
#make count data specific data frame with just the columns you're interested in
count_data <- long.data.2 %>%
dplyr :: select(1:2, 11:12) %>%
  filter(Order!= 'Thysanoptera') #Thysanoptera throws off the data since it is so numerous

##make model with average total insects caught per trap by BurnType
hist(count_data$Count) #data looks negative-binom or Poisson or zero-inflated

##make models and test which is best using AICs and output simulation
#poisson, not zero inflated
count_pois <- glmmTMB(Count ~ BurnType + (1|Forest), data=count_data, 
                      family = "poisson")
#poisson, yes zero inflated
count_pois_zi <- glmmTMB(Count~BurnType + (1|Forest), data=count_data,
                         family="poisson", ziformula=~1)

#negative binomial 1, not zero inflated
count_negbinom1 <- glmmTMB(Count ~ BurnType + (1|Forest), data=count_data,
                           family="nbinom1")

#negative binomial 1, yes zero inflated
count_negbinom1_zi <- glmmTMB(Count ~ BurnType + (1|Forest), data=count_data,
                              family="nbinom1", ziformula=~1)

#negative binomial 2, not zero inflated
count_negbinom2 <- glmmTMB(Count ~ BurnType + (1|Forest), data=count_data,
                           family="nbinom2")

#negative binomial 2, yes zero inflated
count_negbinom2_zi <- glmmTMB(Count ~ BurnType + (1|Forest), data=count_data,
                              family="nbinom2", ziformula=~1)

library(AICcmodavg)
list.1 <- list(count_pois = count_pois, count_pois_zi = count_pois_zi, 
            count_negbinom1 = count_negbinom1, count_negbinom1_zi= count_negbinom1_zi,
            count_negbinom2 = count_negbinom2, count_negbinom2_zi = count_negbinom2_zi)
aictab(list.1)
# Model selection based on AICc:
#   K     AICc Delta_AICc AICcWt Cum.Wt        LL
# count_negbinom2    5  5231.73       0.00   0.54   0.54  -2610.82
# count_negbinom2_zi 6  5233.77       2.04   0.20   0.74  -2610.82
# count_negbinom1    5  5233.84       2.11   0.19   0.93  -2611.88
# count_negbinom1_zi 6  5235.88       4.15   0.07   1.00  -2611.88
# count_pois         4 28351.40   23119.67   0.00   1.00 -14171.67
# count_pois_zi      5 28353.43   23121.70   0.00   1.00 -14171.67
   #count_negbinom2 has the lowest AICc, so that is what we'll go with
simulation_count <- simulateResiduals(fittedModel = count_negbinom2)
plot(simulation_count) #not good. the model does not fit our assumptions for variance and dispersion
testZeroInflation(simulation_count) #p-value is at 0, which means our data is indeed zero inflated

#testing residuals with second lowest AICc, count_negbinom2_zi
simulation_count.1 <- simulateResiduals(fittedModel = count_negbinom2_zi)
plot(simulation_count.1) #not good. the model does not fit our assumptions for variance and dispersion
testZeroInflation(simulation_count.1) # p-value is still 0

#testing residuals for third lowest AICc, count_negbinom1
simulation_count.2 <- simulateResiduals(fittedModel = count_negbinom1)
plot(simulation_count.2) #not good. the model does not fit our assumptions for variance and dispersion
testZeroInflation(simulation_count.2) #p-value still 0

#testing residuals for fourth lowest AICc, count_negbinom1_zi
simulation_count.3 <- simulateResiduals(fittedModel = count_negbinom1_zi)
plot(simulation_count.3) #not good. the model does not fit our assumptions for variance and dispersion
testZeroInflation(simulation_count.3) #p-value still 0

#testing residuals for fifth lowest AICc, count_pois 
simulation_count.4 <- simulateResiduals(fittedModel = count_pois)
plot(simulation_count.4) #not good. the model does not fit our assumptions for variance and dispersion
testZeroInflation(simulation_count.4) #p-value is 1. 

#just going to test the last model, even though it probably will not work
simulation_count.5 <- simulateResiduals(fittedModel = count_pois_zi)
plot(simulation_count.5) #not good. the model does not fit our assumptions for variance and dispersion
testZeroInflation(simulation_count.5) #p-value is 1

####Figures of Data####
##stacked bar graph showing orders and their abundances within forest-sites
long.data.2$BurnType <- factor(long.data.2$BurnType, 
                         levels = c("recent", "medium", "unburned"))

#make new column with combined Forest and BurnType
data_orders <- long.data.2 %>%
  filter(Order != 'Thysanoptera') %>%
  mutate(forestsite = paste(Forest, BurnType, sep = "-"))

all_orders <- ggplot(data_orders, aes(x=forestsite, y = Count)) + 
  geom_bar(stat="identity", aes(fill = Order)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
all_orders

##stacked bar graph showing all orders grouped by burn type alone
burn_orders <- ggplot(data_orders, aes(x=BurnType, y = Count)) + 
  geom_bar(stat="identity", aes(fill = Order)) + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))
burn_orders


###Hymenoptera####
hym_data <- count_data %>%
  filter(Order == "Hymenoptera")

hist(hym_data$Count)

hym_model.gaus<-glmmTMB(Count ~ BurnType + (1|Forest), data = hym_data, 
                        family = "gaussian")
hym_model.nb1 <- update(hym_model.gaus, family = "nbinom1")
hym_model.nb2 <- update(hym_model.gaus, family = "nbinom2")

hym_model.nb2 <- glmmTMB(Count ~ BurnType + (1|Forest), data = hym_data,
                         family = "nbinom2")
hym_model.nb1.zi <- update(hym_model.nb1, ziformula = ~1) 
# hym_model.nb2.zi <- glmmTMB(Count ~ BurnType + (1|Forest), data = hym_data,
#                             family = "nbinom2", ziformula = ~1)
hym_model.pois <- update(hym_model.gaus, family = "poisson") 
hym_model.pois.zi <-  update(hym_model.gaus, family = "poisson", ziformula=~1)

list <- list(hym_model.gaus = hym_model.gaus, hym_model.nb1=hym_model.nb1,
             hym_model.nb2=hym_model.nb2, hym_model.nb1.zi=hym_model.nb1.zi,
             hym_model.pois=hym_model.pois, hym_model.pois.zi=hym_model.pois.zi)
aictab(list)#lowest aic was for hym_model.nb2
hym_model <- hym_model.nb2
summary(hym_model) #see the AIC, not of the null model though

#test for normality metrics. All outputs should be non-significant
simulation_hym <- simulateResiduals(fittedModel = hym_model)
plot(simulation_hym) #meets out assumptions; we are good to continue

 #run the ANOVA on the glm to see the results
hym_anova<-Anova(hym_model, type = "III")
hym_anova

#make summary data containing just the means and standard erros for count data
library(tidyverse)
library(plotrix)
library(ggthemes)
hym_summary <- hym_data %>%
  dplyr::group_by(BurnType) %>%
  dplyr::summarise(
    mean = mean(Count, na.rm = TRUE),
    se = std.error(Count, na.rm = TRUE),
    .groups = "keep"
  )

# Emmeans for treatments
library(emmeans)
library(multcomp)
library(multcompView)

# Get estimated marginal means (emmeans) for the Treatment factor
hym_emmeans <- emmeans(hym_model, specs = pairwise ~ BurnType,
                       adjust = "sidak")
hym_cld_results <- cld(hym_emmeans, adjust = "sidak", Letters = letters)
hym_letters_df <- as.data.frame(hym_cld_results)
print(names(hym_letters_df))

hym_letters_df <- hym_letters_df %>%
  dplyr::select(BurnType, .group) %>%
  dplyr::rename(Letter = .group)

# Merge letters with summary data
hym_summary <- hym_summary %>%
  left_join(hym_letters_df, by = "BurnType")
glimpse(hym_summary)

hym_summary$BurnType <- factor(hym_summary$BurnType, 
                               levels = c("recent", "medium", "unburned"))
hym_graph <- hym_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 6, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE) + 
  ggtitle("Hymenoptera")
hym_graph

hym_data %>%
  mutate(BurnType = factor(BurnType, levels = c("recent", "medium", "unburned")))%>%
  ggplot(., aes(x = BurnType, y = Count, fill = BurnType, 
                color = BurnType, label=Forest))+
  geom_point()+
  geom_text(hjust=0, vjust=0)+
  theme(legend.position="none")+
  ggtitle("Hymenoptera")

hym_graph.poster <- hym_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 7, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE) +
  scale_fill_manual(values = c("#c6dbef", "#6baed6",
                               "#08306b"))+
  scale_x_discrete(labels = c("Recent Burn", "Intermediate Burn", "Unburned"))+
  labs(x = "Prescribed Fire History", y = "Mean Abundance Per Trap")+
  theme_clean()+
  theme(axis.title.x = element_text(size = 15, family = "Tahoma"),
        axis.title.y = element_text(size = 15, family = "Tahoma"),
        axis.text.x = element_text(size = 12, family = "Tahoma"),
        axis.text.y = element_text(size = 12, family = "Tahoma"),
        legend.text = element_text(family = "Tahoma"),
        legend.title = element_text(family = "Tahoma"),
        plot.title = element_text(size = 25, family = "Tahoma",
                                  hjust = 0.85, vjust = -1.3))+
  ggtitle("Hymenoptera")
hym_graph.poster

#Coleoptera####
col_data <- count_data %>%
  filter(Order == "Coleoptera")

hist(col_data$Count)

col_model.gaus<-glmmTMB(Count ~ BurnType + (1|Forest), data = col_data, 
                        family = "gaussian")
col_model.nb1 <- update(col_model.gaus, family = "nbinom1")
col_model.nb2 <- update(col_model.gaus, family = "nbinom2")
col_model.nb1.zi <- update(col_model.nb1, ziformula = ~1) 
col.model.nb2.zi <- update(col_model.nb2, ziformula=~1)
col_model.pois <- update(col_model.gaus, family = "poisson") 
col_model.pois.zi <-  update(col_model.gaus, family = "poisson", ziformula=~1)

list <- list(col_model.gaus = col_model.gaus, col_model.nb1=col_model.nb1,
             col_model.nb2=col_model.nb2, col_model.nb1.zi=col_model.nb1.zi,
             col_model.pois=col_model.pois, col_model.pois.zi=col_model.pois.zi,
             col.model.nb2.zi = col.model.nb2.zi)
aictab(list) #lowest aic was for col_model.nb1, but the assumptions fit better for col_model.nb2
col_model <- col_model.nb2
summary(col_model) 

#test for normality metrics. All outputs should be non-significant
simulation_col <- simulateResiduals(fittedModel = col_model)
plot(simulation_col) #some deviation. does not fit out assumptions so well



#run the ANOVA on the glm to see the results
col_anova<-Anova(col_model, type = "III")
col_anova

#make summary data containing just the means and standard erros for count data
col_summary <- col_data %>%
  dplyr::group_by(BurnType) %>%
  dplyr::summarise(
    mean = mean(Count, na.rm = TRUE),
    se = std.error(Count, na.rm = TRUE),
    .groups = "keep"
  )

# Get estimated marginal means (emmeans) for the Treatment factor
col_emmeans <- emmeans(col_model, specs = pairwise ~ BurnType,
                       adjust = "sidak")
col_cld_results <- cld(col_emmeans, adjust = "sidak", Letters = letters)
col_letters_df <- as.data.frame(col_cld_results)
print(names(col_letters_df))

col_letters_df <- col_letters_df %>%
  dplyr::select(BurnType, .group) %>%
  dplyr::rename(Letter = .group)

# Merge letters with summary data
col_summary <- col_summary %>%
  left_join(col_letters_df, by = "BurnType")
glimpse(col_summary)

col_summary$BurnType <- factor(col_summary$BurnType, 
                               levels = c("recent", "medium", "unburned"))
col_graph <- col_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 6, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE) +
  ggtitle("Coleoptera")
col_graph

col_data %>%
  mutate(BurnType = factor(BurnType, 
                           levels = c("recent", "medium", "unburned")))%>%
  ggplot(., aes(x = BurnType, y = Count, fill = BurnType, 
                color = BurnType, label = Forest))+
  geom_point()+
  geom_text(hjust=0, vjust=0)+
  theme(legend.position="none")+
  ggtitle("Coleoptera")

col_graph.poster <- col_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 7, color = "black",
            hjust = 0.5, vjust = -0.1,
            position = position_nudge(x=0),
            parse = TRUE) +
  scale_fill_manual(values = c("#fc9272", "#ef3b2c",
                               "#99000d"))+
  scale_x_discrete(labels = c("Recent Burn", "Intermediate Burn", "Unburned"))+
  labs(x = "Prescribed Fire History", y = "Mean Abundance Per Trap")+
  theme_clean()+
  theme(axis.title.x = element_text(size = 15, family = "Tahoma"),
        axis.title.y = element_text(size = 15, family = "Tahoma"),
        axis.text.x = element_text(size = 12, family = "Tahoma"),
        axis.text.y = element_text(size = 12, family = "Tahoma"),
        legend.text = element_text(family = "Tahoma"),
        legend.title = element_text(family = "Tahoma"),
        plot.title = element_text(size = 25, family = "Tahoma",
                                  hjust = 0.85, vjust = -1.3))+
  ggtitle("Coleoptera")
col_graph.poster

#Hemiptera####
hem_data <- count_data %>%
  filter(Order == "Hemiptera")

hist(hem_data$Count)

hem_model.gaus<-glmmTMB(Count ~ BurnType + (1|Forest), data = hem_data, 
                        family = "gaussian")

hem_model.nb1 <- update(hem_model.gaus, family = "nbinom1")
hem_model.nb2 <- update(hem_model.gaus, family = "nbinom2")
hem_model.nb1.zi <- update(hem_model.nb1, ziformula = ~1) 
hem_model.nb2.zi <- update(hem_model.nb2, ziformula=~1)
hem_model.pois <- update(hem_model.gaus, family = "poisson") 
hem_model.pois.zi <-  update(hem_model.gaus, family = "poisson", 
                             ziformula=~1)



list <- list(hem_model.gaus = hem_model.gaus, hem_model.nb1=hem_model.nb1,
             hem_model.nb2=hem_model.nb2, hem_model.nb1.zi=hem_model.nb1.zi,
             hem_model.pois=hem_model.pois, hem_model.pois.zi=hem_model.pois.zi,
             hem_model.nb2.zi = hem_model.nb2.zi)
aictab(list) #lowest aic was for hem_model.nb1
hem_model <- hem_model.nb1

summary(hem_model) #see the AIC, not of the null model though

#test for normality metrics. All outputs should be non-significant
simulation_hem <- simulateResiduals(fittedModel = hem_model)
plot(simulation_hem) #some within group deviation. otherwise assumptions are met

#run the ANOVA on the glm to see the results
hem_anova<-Anova(hem_model, type = "III")
hem_anova

#make summary data containing just the means and standard erros for count data
hem_summary <- hem_data %>%
  dplyr::group_by(BurnType) %>%
  dplyr::summarise(
    mean = mean(Count, na.rm = TRUE),
    se = std.error(Count, na.rm = TRUE),
    .groups = "keep"
  )

# Get estimated marginal means (emmeans) for the Treatment factor
hem_emmeans <- emmeans(hem_model, specs = pairwise ~ BurnType,
                       adjust = "sidak")
hem_cld_results <- cld(hem_emmeans, adjust = "sidak", Letters = letters)
hem_letters_df <- as.data.frame(hem_cld_results)
print(names(hem_letters_df))

hem_letters_df <- hem_letters_df %>%
  dplyr::select(BurnType, .group) %>%
  dplyr::rename(Letter = .group)

# Merge letters with summary data
hem_summary <- hem_summary %>%
  left_join(hem_letters_df, by = "BurnType")
glimpse(hem_summary)

hem_summary$BurnType <- factor(hem_summary$BurnType, 
                               levels = c("recent", "medium", "unburned"))
hem_graph <- hem_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 6, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE) +
  ggtitle("Hemiptera")
hem_graph

hem_data %>%
  mutate(BurnType = factor(BurnType, levels = c("recent", "medium", "unburned")))%>%
  ggplot(., aes(x = BurnType, y = Count, 
                fill = BurnType, color = BurnType, label = Forest))+
  geom_point()+
  geom_text(hjust=0, vjust=0)+
  theme(legend.position="none")+
  ggtitle("Hemiptera")


hem_graph.poster <- hem_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 7, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE) +
  scale_fill_manual(values = c("#dadaeb", "#807dba",
                               "#3f007d"))+
  scale_x_discrete(labels = c("Recent Burn", "Intermediate Burn", "Unburned"))+
  labs(x = "Prescribed Fire History", y = "Mean Abundance Per Trap")+
  theme_clean()+
  ggtitle("Hemiptera")+
  theme(axis.title.x = element_text(size = 15, family = "Tahoma"),
        axis.title.y = element_text(size = 15, family = "Tahoma"),
        axis.text.x = element_text(size = 12, family = "Tahoma"),
        axis.text.y = element_text(size = 12, family = "Tahoma"),
        legend.text = element_text(family = "Tahoma"),
        legend.title = element_text(family = "Tahoma"),
        plot.title = element_text(size = 25, family = "Tahoma",
                                  hjust = 0.85, vjust = -1.3))


hem_graph.poster
#Diptera####
dip_data <- count_data %>%
  filter(Order == "Diptera")

hist(dip_data$Count)

dip_model.gaus<-glmmTMB(Count ~ BurnType + (1|Forest), data = dip_data, 
                        family = "gaussian")

dip_model.nb1 <- update(dip_model.gaus, family = "nbinom1")
dip_model.nb2 <- update(dip_model.gaus, family = "nbinom2")
dip_model.nb1.zi <- update(dip_model.nb1, ziformula = ~1) 
dip_model.nb2.zi <- update(dip_model.nb2, ziformula=~1)
dip_model.pois <- update(dip_model.gaus, family = "poisson") 
dip_model.pois.zi <-  update(dip_model.gaus, family = "poisson", 
                             ziformula=~1)



list <- list(dip_model.gaus = dip_model.gaus, dip_model.nb1=dip_model.nb1,
             dip_model.nb2=dip_model.nb2, dip_model.nb1.zi=dip_model.nb1.zi,
             dip_model.pois=dip_model.pois, dip_model.pois.zi=dip_model.pois.zi,
             dip_model.nb2.zi = dip_model.nb2.zi)
aictab(list) #lowest aic was for dip_model.nb1
dip_model <- dip_model.nb1
summary(dip_model) #see the AIC, not of the null model though

#test for normality metrics. All outputs should be non-significant
simulation_dip <- simulateResiduals(fittedModel = dip_model)
plot(simulation_dip) #NICE! All assumptions are met

#run the ANOVA on the glm to see the results
dip_anova<-Anova(dip_model, type = "III")
dip_anova

#make summary data containing just the means and standard erros for count data
dip_summary <- dip_data %>%
  dplyr::group_by(BurnType) %>%
  dplyr::summarise(
    mean = mean(Count, na.rm = TRUE),
    se = std.error(Count, na.rm = TRUE),
    .groups = "keep"
  )

# Get estimated marginal means (emmeans) for the Treatment factor
dip_emmeans <- emmeans(dip_model, specs = pairwise ~ BurnType,
                       adjust = "sidak")
dip_cld_results <- cld(dip_emmeans, adjust = "sidak", Letters = letters)
dip_letters_df <- as.data.frame(dip_cld_results)
print(names(dip_letters_df))

dip_letters_df <- dip_letters_df %>%
  dplyr::select(BurnType, .group) %>%
  dplyr::rename(Letter = .group)

# Merge letters with summary data
dip_summary <- dip_summary %>%
  left_join(dip_letters_df, by = "BurnType")
glimpse(dip_summary)

dip_summary$BurnType <- factor(dip_summary$BurnType, 
                               levels = c("recent", "medium", "unburned"))
dip_graph <- dip_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 6, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE)+
  ggtitle ("Diptera")
dip_graph

dip_data %>%
  mutate(BurnType = factor(BurnType, levels = c("recent", "medium", "unburned")))%>%
  ggplot(., aes(x = BurnType, y = Count, fill = BurnType, 
                color = BurnType, label = Forest))+
  geom_point()+
  geom_text(hjust=0,vjust=0)+
  theme(legend.position="none")+
  ggtitle("Diptera")

dip_graph.poster <- dip_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 7, color = "black",
            hjust = 0.5, vjust = -0.2,
            position = position_nudge(x=0),
            parse = TRUE) +
  scale_fill_manual(values = c("#fec44f", "#ec7014",
                               "#993404"))+
  scale_x_discrete(labels = c("Recent Burn", "Intermediate Burn", "Unburned"))+
  labs(x = "Prescribed Fire History", y = "Mean Abundance Per Trap")+
  theme_clean()+
  theme(axis.title.x = element_text(size = 15, family = "Tahoma"),
        axis.title.y = element_text(size = 15, family = "Tahoma"),
        axis.text.x = element_text(size = 12, family = "Tahoma"),
        axis.text.y = element_text(size = 12, family = "Tahoma"),
        legend.text = element_text(family = "Tahoma"),
        legend.title = element_text(family = "Tahoma"))+
  ggtitle("Diptera")

dip_graph.poster

###Richness####
#use the data frame with each order as its own column
  #but first, remove the Thysanoptera since you stopped collecting this data
wide.data.2 <- wide.data.1 %>%
  dplyr::select(1:23, 25) #some other package in use has a select() function

data.rich <- wide.data.2 %>%
  mutate(Richness = apply(wide.data.2[,7:23]>0, 1, sum))

#visualize the data
hist(data.rich$Richness)
  #looks like a normal (Gaussian) distribution
rich_model.gaus<-glmmTMB(Richness ~ BurnType + (1|Forest), data = data.rich, 
                         family = "gaussian")

rich_model.nb1 <- update(rich_model.gaus, family = "nbinom1")
rich_model.nb2 <- update(rich_model.gaus, family = "nbinom2")
rich_model.nb1.zi <- update(rich_model.nb1, ziformula = ~1) 
rich_model.nb2.zi <- update(rich_model.nb2, ziformula=~1)
rich_model.pois <- update(rich_model.gaus, family = "poisson") 
rich_model.pois.zi <-  update(rich_model.gaus, family = "poisson", 
                              ziformula=~1)
rich_model.tweed <- update(rich_model.gaus, family = "tweedie")



list <- list(rich_model.gaus = rich_model.gaus, rich_model.nb1=rich_model.nb1,
             rich_model.nb2=rich_model.nb2, rich_model.nb1.zi=rich_model.nb1.zi,
             rich_model.pois=rich_model.pois, rich_model.pois.zi=rich_model.pois.zi,
             rich_model.nb2.zi = rich_model.nb2.zi,
             rich_model.tweed = rich_model.tweed)
aictab(list) #lowest AIC is rich_model.gaus. Probably because we don't have Count data
#and count data is used for the other distribution types. Also, the distribution

rich_model <- rich_model.gaus

summary(rich_model) #see the AIC, not of the null model though

#test for normality metrics. All outputs should be non-significant
simulation_rich <- simulateResiduals(fittedModel = rich_model)
plot(simulation_rich) #nothing is p<0.05, which means all assumptions are met

#run the ANOVA on the glm to see the results
rich_anova<-Anova(rich_model, type = "III")
rich_anova #p-value of 0.2133, appears to be no significance

#make summary data containing just the means and standard errors for count data
rich_summary <- data.rich %>%
  dplyr::group_by(BurnType) %>%
  dplyr::summarise(
    mean = mean(Richness, na.rm = TRUE),
    se = std.error(Richness, na.rm = TRUE),
    .groups = "keep"
  )

# Get estimated marginal means (emmeans) for the Treatment factor
rich_emmeans <- emmeans(rich_model, specs = pairwise ~ BurnType,
                        adjust = "sidak")
rich_cld_results <- cld(rich_emmeans, adjust = "sidak", Letters = letters)
rich_letters_df <- as.data.frame(rich_cld_results)
print(names(rich_letters_df))

rich_letters_df <- rich_letters_df %>%
  dplyr::select(BurnType, .group) %>%
  dplyr::rename(Letter = .group)

# Merge letters with summary data
rich_summary <- rich_summary %>%
  left_join(rich_letters_df, by = "BurnType")
glimpse(rich_summary)

rich_summary$BurnType <- factor(rich_summary$BurnType, 
                                levels = c("recent", "medium", "unburned"))
rich_graph <- rich_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 6, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE) +
  ggtitle("Richness")
rich_graph #finding: no significant differences in richness

#one more graph so you can see individual data points, observe effect of forest
data.rich %>%
  mutate(BurnType = factor(BurnType, levels = c("recent", "medium", "unburned")))%>%
  ggplot(., aes(x = BurnType, y = Richness, 
                fill = BurnType, color = BurnType, label = Forest))+
  geom_point()+
  geom_text(hjust=0, vjust=0)+
  theme(legend.position="none")+
  ggtitle("Richness")

##Evenness####
#Evenness Formula:
#Pielou's Evenness = Shannon Diversity Index / ln(richness)
data.ev <- data.rich %>%
  mutate(shannon = data_shan$shan) %>%
  mutate(Evenness = data_shan$shan / log(data.rich$Richness))

#visualize the data
hist(data.ev$Evenness) #looks normal (Gaussian) distribution

ev_model.gaus<-glmmTMB(Evenness ~ BurnType + (1|Forest), data = data.ev, 
                       family = "gaussian")

ev_model.nb1 <- update(ev_model.gaus, family = "nbinom1")
ev_model.nb2 <- update(ev_model.gaus, family = "nbinom2")
ev_model.nb1.zi <- update(ev_model.nb1, ziformula = ~1) 
ev_model.nb2.zi <- update(ev_model.nb2, ziformula=~1)
ev_model.pois <- update(ev_model.gaus, family = "poisson") 
ev_model.pois.zi <-  update(ev_model.gaus, family = "poisson", 
                            ziformula=~1)
ev_model.tweed <- update(ev_model.gaus, family = "tweedie")



list <- list(ev_model.gaus = ev_model.gaus, ev_model.nb1=ev_model.nb1,
             ev_model.nb2=ev_model.nb2, ev_model.nb1.zi=ev_model.nb1.zi,
             ev_model.pois=ev_model.pois, ev_model.pois.zi=ev_model.pois.zi,
             ev_model.nb2.zi = ev_model.nb2.zi,
             ev_model.tweed = ev_model.tweed)
aictab(list) #lowest AIC is Guassian

ev_model <- ev_model.gaus

summary(ev_model) #see the AIC, not of the null model though

#test for normality metrics. All outputs should be non-significant
simulation_ev <- simulateResiduals(fittedModel = ev_model)
plot(simulation_ev) #no p-values below 0.05, all assumptions met

#run the ANOVA on the glm to see the results
ev_anova<-Anova(ev_model, type = "III")
ev_anova #p-value of 0.8238, very non-significant

#make summary data containing just the means and standard erros for count data
ev_summary <- data.ev %>%
  dplyr::group_by(BurnType) %>%
  dplyr::summarise(
    mean = mean(Evenness, na.rm = TRUE),
    se = std.error(Evenness, na.rm = TRUE),
    .groups = "keep"
  )

# Get estimated marginal means (emmeans) for the Treatment factor
ev_emmeans <- emmeans(ev_model, specs = pairwise ~ BurnType,
                      adjust = "sidak")
ev_cld_results <- cld(ev_emmeans, adjust = "sidak", Letters = letters)
ev_letters_df <- as.data.frame(ev_cld_results)
print(names(ev_letters_df))

ev_letters_df <- ev_letters_df %>%
  dplyr::select(BurnType, .group) %>%
  dplyr::rename(Letter = .group)

# Merge letters with summary data
ev_summary <- ev_summary %>%
  left_join(ev_letters_df, by = "BurnType")
glimpse(ev_summary)

ev_summary$BurnType <- factor(ev_summary$BurnType, 
                              levels = c("recent", "medium", "unburned"))
ev_graph <- ev_summary %>%
  group_by(BurnType) %>%
  ggplot(., aes(x = BurnType, y = mean, fill = BurnType, color = BurnType))+
  geom_col(show.legend = FALSE, color = "black", alpha = 0.8)+
  geom_errorbar(aes(ymin = mean - se, ymax = mean+se),
                width = 0.05, size = 0.5, color = "black")+
  geom_text(aes(label = Letter, 
                y = mean + se + 0.5),
            size = 6, color = "black",
            hjust = 0.5, vjust = 0,
            position = position_nudge(x=0),
            parse = TRUE) +
  ggtitle("Evenness")
ev_graph #finding: no significant differences in Evenness

#another graph to see individual points
data.ev %>%
  mutate(BurnType = factor(BurnType, levels = c("recent", "medium", "unburned")))%>%
  ggplot(., aes(x = BurnType, y = Evenness, 
                fill = BurnType, color = BurnType, label = Forest))+
  geom_point()+
  geom_text(hjust=0, vjust=0)+
  theme(legend.position="none")+
  ggtitle("evness")

###NMDS Plot####
library(vegan)
#This code is originally from Danilo Dos Santos based on his data from a project from Summer 2023


# First, pull out the columns in your dataframe that are just your counts. 
#Your resulting dataset should only have the Orders and their counts.
community <- wide.data.2[,7:23]

# This runs the Hellinger transformation on the data, which accounts for 
#"rare" species. It's a form of scaling.
comm.hel <- decostand(community,"hellinger")


# This turns that into a matrix, which is needed for the NMDS.
m_com = as.matrix(comm.hel)

# this runs the NMDS. Since I don't specify a distance, it's running the 
#Bray-Curtis, but if you think there's a distance that would work better for 
#your data, you can use another. Full list: "manhattan", "euclidean", "canberra", 
#"clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", 
#"horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq",
#"chord", "hellinger", "aitchison", or "robust.aitchison"
nmds = metaMDS(m_com) 

# Shows the NMDS data. You really want to look at your stress here. You want the stress to be below 0.2. 
nmds #stress is 0.234, which is not ideal. But it is nearly below 0.2. Not sure what to do with this.
stressplot(nmds) #points roughly follow along the stress plot

# This creates an NMDS plot. First line creates your blank plot and puts in 
#what you specify in type. Type = "t" for text, "p" for points, or "n" for none. 
#Second line displays your species (or Family or Order) data, with red text and making sure 
#there's space between labels (air). Third line displays your sites with black 
#text (default) with text size 1.25 and air.
ordiplot(nmds, type="n")
orditorp(nmds, display="species",col="red",air=0.01)
orditorp(nmds, display="sites",cex=1.25,air=0.01)

# This creates the fun graph where you can overlay polygons for burn type or site!
#First line I'm calling a new plot and I like point so I choose points. 
#Ordihull would give you a polygon that has hard edges - play around with it 
#to see if you like it. I've chosen specific colors for my 3 unique trap types, 
#but you can put your own colors in. They go in the order that they were in 
#your original dataset (comm.july2023). I also like a black border on my 
#polygons. Orditorp is just bringing in my species names so I can see where 
#they are in orthogonal space. Ordiellipse is the circles around your data 
#point - for me here it's Trap Type. You could also call comm.july2023$Site 
#instead (for you it'd be Forest) to see if your communities were distinct across your forests. 

ordiplot(nmds,type="p")
ordihull(nmds,groups=wide.data.2$BurnType,draw="polygon",
         col=c("#e9de4f", "#2e6bb6", "#aaa7a7"),label=TRUE, border = "black")
orditorp(nmds,display="species",col="black",air=0.01)
ordiellipse(nmds,groups=wide.data.2$BurnType,draw="polygon", col=c("#e9de4f", "#2e6bb6", "#aaa7a7"), label=FALSE,border="black")
  #looks like there is a whole lot of overlap here. I would say no significance.
# Adonis2 is the standard vegan() for assessing if your community differs by 
#your chosen variable (for me, TrapType). Check the p-value, simple as that.
adonis2(comm.hel~BurnType, data=wide.data.2)

# I (Danilo) found a package online that does a pairwise adonis2() on the data to check 
#the differences - there might be a better way to do this idk? If you want to 
#download it, you have to get it from the guy's github, but you'd also have to 
#clear your github token on R blahblahblah. it was a pain. I would see if 
#anyone else has a better way to do it. 

#pairwise.adonis2(comm.hel~TrapType, data=comm.july2023, nperm=999)

# Here I run an indicator species analysis (package indicspecies, function multipatt).
#It allows determining lists of species that are associated to particular 
#groups of sites (or combinations of those). 

library(indicspecies)

# First, redo the decostand() on your basic community data, but this time do presence/absence.
pa <- decostand(community, "pa")

# this then creates the indicator analysis and tells it what to cluster by.
burn.pa <- multipatt(x = pa, cluster = wide.data.2$BurnType)

# look at that stuff. If there is something that is significantly indicated for 
#a specific trap, it will show in its own little anova with a p-value. 
#If there's nothing, then nothing was indicated to be associated with a particular trap.
summary(burn.pa)
