library(dplyr)
library(tidyr)
library(rstatix)
library(ggplot2)


load("s_u.rdata") #urine
load("s_p.rdata") #plasma
load("s_k.rdata") #kidney
load("s_h.rdata") #heart
load("s_l.rdata") #liver

#process and statistic analysis
#s_u
s_u_e <- s_u %>% #data
  group_by(Compound) %>%
  mutate(n = n()) %>%
  group_by(Compound, group) %>%
  mutate(n2 = n()) %>%
  filter(n2 <= 1)

s_u_ev <- unique(s_u_e$Compound) 

s_u_f <- s_u %>% #data
  filter(!Compound %in% c(s_u_ev, c("Oxaloacetate", "Inosine", "G6P"))) %>% #Remove non-lysine/TCA metabolites or single-sample metabolites
  mutate(log2_Area = log2(Area)) %>%
  group_by(Compound) %>%
  mutate(group = factor(group, levels = c("KO", "WT"))) %>%
  wilcox_test(log2_Area~group, detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  select(Compound, fc = estimate, p, p.adj) %>%
  mutate(Model = "Mlycd KO", #model
         Organ = "Urine")#organ

remove(s_u_e, s_u_ev)

#s_p
s_p_e <- s_p %>% #data
  group_by(Compound) %>%
  mutate(n = n()) %>%
  group_by(Compound, group) %>%
  mutate(n2 = n()) %>%
  filter(n <= 4)

s_p_ev <- unique(s_p_e$Compound) 

s_p_f <- s_p %>% #data
  filter(!Compound %in% c(s_p_ev, c("Glutamic acid", "Inosine", "G6P"))) %>% #Remove non-lysine/TCA metabolites or single-sample metabolites
  mutate(log2_Area = log2(Area)) %>%
  group_by(Compound) %>%
  mutate(group = factor(group, levels = c("KO", "WT"))) %>%
  wilcox_test(log2_Area~group, detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  select(Compound, fc = estimate, p, p.adj) %>%
  mutate(Model = "Mlycd KO", #model
         Organ = "Plasma")#organ

remove(s_p_e, s_p_ev)

#s_k
s_k_e <- s_k %>% #data
  group_by(Compound) %>%
  mutate(n = n()) %>%
  group_by(Compound, group) %>%
  mutate(n2 = n()) %>%
  filter(n2 <=1) 

s_k_ev <- unique(s_k_e$Compound) 

s_k_f <- s_k %>% #data
  filter(!Compound %in% c(s_k_ev, c("Glutamic acid", "Inosine", "G6P", "GAP"))) %>% #Remove non-lysine/TCA metabolites or single-sample metabolites
  mutate(log2_Area = log2(Area)) %>%
  group_by(Compound) %>%
  mutate(group = factor(group, levels = c("KO", "WT"))) %>%
  wilcox_test(log2_Area~group, detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  select(Compound, fc = estimate, p, p.adj) %>%
  mutate(Model = "Mlycd KO", #model
         Organ = "Kidney")#organ

remove(s_k_e, s_k_ev)

#s_l
s_l_e <- s_l %>% #data
  group_by(Compound) %>%
  mutate(n = n()) %>%
  group_by(Compound, group) %>%
  mutate(n2 = n()) %>%
  filter(n2 <=1 ) 

s_l_ev <- unique(s_l_e$Compound) 

s_l_f <- s_l %>% #data
  filter(!Compound %in% c(s_l_ev, c("Glutamic acid", "Inosine", "G6P", "GAP", "Glutathione", "F16P", "2PG", "3PG"))) %>% #Remove non-lysine/TCA metabolites or single-sample metabolites
  mutate(log2_Area = log2(Area)) %>%
  group_by(Compound) %>%
  mutate(group = factor(group, levels = c("KO", "WT"))) %>%
  wilcox_test(log2_Area~group, detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  select(Compound, fc = estimate, p, p.adj) %>%
  mutate(Model = "Mlycd KO", #model
         Organ = "Liver")  #organ 

remove(s_l_e, s_l_ev)

#heart
s_h_e <- s_h %>% #data
  group_by(Compound) %>%
  mutate(n = n()) %>%
  group_by(Compound, group) %>%
  mutate(n2 = n()) %>% 
  filter(n <= 1)

s_h_ev <- unique(s_h_e$Compound) 

s_h_f <- s_h %>% #data
  filter(!Compound %in% c(s_h_ev, c("Glutamic acid", "Inosine", "G6P", "3PG", "Glutathione", "GAP"))) %>% #Remove non-lysine/TCA metabolites or single-sample metabolites
  mutate(log2_Area = log2(Area)) %>%
  group_by(Compound) %>%
  mutate(group = factor(group, levels = c("KO", "WT"))) %>%
  wilcox_test(log2_Area~group, detailed = T) %>%
  adjust_pvalue(method = "fdr") %>%
  select(Compound, fc = estimate, p, p.adj) %>%
  mutate(Model = "Mlycd KO", #model
         Organ = "Heart")#organ

remove(s_h_e, s_h_ev)


data <- rbind(s_k_f, s_l_f, s_p_f, s_u_f, s_h_f) %>%
  mutate(Organ = factor(Organ, levels = c("Heart", "Kidney", "Liver", "Plasma", "Urine")),
         sig = case_when(
           p < 0.001 ~ "***",
           p < 0.01 ~ "**",
           p < 0.05 ~ "*",
           T ~ ""
         ),
         Compound = case_when(
           Compound == "2PG" | Compound == "3PG" ~ "2PG/3PG",
           Compound == "F6P" ~ "G6P/F6P",
           Compound == "DHAP" ~ "GAP/DHAP",
           Compound == "alpha-Ketoglutarate" ~ "α-Ketoglutarate",
           T ~ Compound
         ),
         Compound = factor(Compound, levels = c("Glucose", "G6P/F6P","F16P", "GAP/DHAP", "2PG/3PG",
                                                "PEP", "Lactate", "Pyruvate", "Citrate", "cis-Aconitate", "Isocitrate","α-Ketoglutarate",
                                                "Succinate","Fumarate","Malate")),
         Compound = fct_rev(Compound)) %>%
  complete(Compound, Organ)

p <- ggplot(data, aes(x = Organ, y = Compound, fill = fc)) +
  
  geom_tile(color = "black", linewidth = 0.5) +
  geom_text(aes(label = sig), color = "black", size = 8, vjust = 0.75) +
  #geom_raster()+
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),      
    values = scales::rescale(c(-5, 0, 2)), 
    limits = c(-5, 2),
    na.value = "grey80"
  ) +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color = "black", size = 24),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 24, color = "black"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)
  ) +
  coord_fixed(ratio = 1) +
  labs(fill = "Log2FC") 