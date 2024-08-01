library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(rlang)
library(car)

method2_anova <- function(df, x){
  
  # get data the same length
  pre_df = df %>%
    select({{x}}, Feedback, technique, inflation, Participant)%>%
    group_by(Feedback) 
  
  counts =  pre_df %>% group_by(Feedback, technique) %>% 
    count(Feedback) %>% 
    summarise(min(n))
  n_1_tech = min(counts[counts$technique == 'single_person',]['min(n)'])
  n_2_tech = min(counts[counts$technique == 'two_person',]['min(n)'])
  # ( n_1_tech, n_2_tech)
  pair_df = map2_df(unique(pre_df$technique), c( n_1_tech, n_2_tech),
                    ~pre_df %>% 
                      filter(technique == .x) %>% 
                      slice_sample(n = .y))
  
  # Summary results
  summary_stats <- pair_df%>%
    group_by(Feedback, technique) %>%
    get_summary_stats(x, type = "full")

  
  #Shapiro-Wilk normality test
  sw <- pair_df %>%
    group_by(technique, Feedback) %>%
    shapiro_test(!!sym(x))



  
  # Pairwise t-test
  pwc1 <- pair_df %>%
    select(condition_var = {{x}}, technique, Feedback)%>%
    group_by(technique) %>%
    pairwise_t_test(
      condition_var ~ Feedback, paired = TRUE, detailed = TRUE,
      p.adjust.method = "bonferroni"
    )

  # Wilcoxon test
  pwc <- pair_df  %>% 
    select(condition_var = {{x}}, technique, Feedback)%>%
    group_by(technique) %>%
    wilcox_test(condition_var ~Feedback, paired = TRUE,detailed = TRUE,) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()
  
  # ANOVA
  wdat <- tidyr::pivot_wider(pair_df,
                             id_cols = c("Participant", "inflation", "technique"),
                             names_from = Feedback, values_from = x, values_fn = mean)
  
  
  wdat$technique  <- as.factor(wdat$technique )
  # a linear model for a multivariate dependent variable
  mvmod <- lm(cbind(F, NF) ~ technique, data=wdat)
  
  # construct an object which reflects the structure
  idata <- data.frame(Feedback = factor(c("F", "NF")))
  # ANOVA
  rmaov <- Anova(mvmod, idata=idata, idesign = ~Feedback,  type=3)
  res = summary(rmaov, multivariate=FALSE)
  aov.res = as.data.frame(res[["univariate.tests"]][,])
  
  
  # Visualisation
  pair_df %>%
    group_by(Feedback, technique) %>%
    get_summary_stats(x, type = "mean_sd")
  
  
  bxp <- ggboxplot(
    pair_df, x = "technique", y = {{x}},
    color = "Feedback", palette = "jco",
    add = "point"
  )
  bxp+theme(legend.position = "right",
            plot.margin = margin(1, 5, 1, 1))
  
  
  # Visualization: box plots with p-values
  # Prepare file path
  # Wilcoxon plot
  # file_name <- sprintf("Plot/ANOVA/ALL/%s_wilcoxon.pdf", x)
  # pwc <- pwc %>% add_xy_position(x = "technique")
  # pdf(file_name,  width = 8, height = 5)
  # print(
  #   bxp +
  #     stat_pvalue_manual(pwc, tip.length = 0, label = "{p} {p.signif}") +
  #     ylab(x) +
  #     labs(
  #       subtitle = (paste("ANOVA p-value:", aov.res["technique:Feedback",6])),
  #       caption = get_pwc_label(pwc)
  #     )
  # )
  # dev.off()
  # # TTest plot
  # file_name <- sprintf("Plot/ANOVA/ALL/%s_ttest.pdf", x)
  # pwc1 <- pwc1 %>% add_xy_position(x = "technique")
  # pdf(file_name,  width = 8, height = 5)
  # print(
  #   bxp +
  #     stat_pvalue_manual(pwc1, tip.length = 0, label = "{p.adj} {p.adj.signif}") +
  #     ylab(x) +
  #     labs(
  #       subtitle = (paste("ANOVA p-value:", aov.res["technique:Feedback",6])),
  #       caption = get_pwc_label(pwc1)
  #     )
  # )
  # dev.off()
  return (list( summary_stats, aov.res, pwc, pwc1, sw))
  # return (list(aov.res, summary_stats, pwc, sw))
}


# Main program
df <- read.csv("Data/processed_data_final.csv")

df$technique[df$technique == 1] <- "single_person"
df$technique[df$technique == 2] <- "two_person"

# leak_diff_per, Inspiratory_time, Expiratory_time, Vti, n(),"Ventilatory_rate"

cols_to_summarize = c("leak_diff_per", "Inspiratory_time", 
                      "Expiratory_time", "Vti")
# Define empty data frame to store output
output_df_aov <- data.frame()
output_df_stat <- data.frame()
output_df_wil <- data.frame()
output_df_pwc <- data.frame()
output_df_sw <- data.frame()

# Main
role_list = c(1,2,3)
mannikin_list = c(1,2)


              
# # ----------------For all the cases



  # Filter to term (2) pr pre-term (1)makinin
mannikin =1
data <- df[df$mannikin == mannikin, ]
  
for (x in cols_to_summarize){
  res = method2_anova(data, x)
  summary_stats = cbind( x, as.data.frame(res[1]))
  wilcox = cbind(mannikin, x, as.data.frame(res[3])[c("technique","estimate",
                                          "group1","group2",
                                          "n1","n2",
                                          "p.adj")])
  pwc = cbind(mannikin, x, as.data.frame(res[4])[c("technique",
                                          "group1","group2",
                                          "n1","n2",
                                          "p.adj")])
  sw = cbind(mannikin, x, as.data.frame(res[5]))
  res_aov = cbind(mannikin, x, as.data.frame(res[2]))
  output_df_aov <- rbind(output_df_aov, res_aov)
  output_df_stat <- rbind(output_df_stat, summary_stats)
  output_df_wil <- rbind(output_df_wil, wilcox)
  output_df_pwc <- rbind(output_df_pwc, pwc)
  output_df_sw <- rbind(output_df_sw, sw)
}



write.csv(output_df_aov, file = "Result/preterm_anova_pvalues.csv")
write.csv(output_df_stat, file = "Result/preterm_anova_summary_stat.csv", row.names = FALSE)
write.csv(output_df_pwc, file = "Result/preterm_tech_ttest_test.csv", row.names = FALSE)
write.csv(output_df_sw, file = "Result/preterm_tech_anova_normality_test.csv", row.names = FALSE)
write.csv(output_df_wil, file = "Result/preterm_tech_wilcoxon_test_test.csv", row.names = FALSE)


# Post-hoc Turkey
# Select the variable
library(emmeans)
library(afex)

# Make it to the same length
pre_df = df %>%
  # select(x, Feedback, technique, inflation, Participant)%>%
  group_by(Feedback)%>%
  arrange(technique)
counts =  pre_df %>% group_by(Feedback, technique) %>% 
  count(Feedback) %>% 
  summarise(min(n))
n_1_tech = min(counts[counts$technique == 'single_person',]['min(n)'])
n_2_tech = min(counts[counts$technique == 'two_person',]['min(n)'])

pair_df = map2_df(unique(pre_df$technique), c(n_1_tech, n_2_tech),
                  ~pre_df %>% 
                    filter(technique == .x) %>% 
                    slice_sample(n = .y)) 

# leak_diff_per
afmod <- afex::aov_car(Vti ~ Feedback * technique 
                       + Error(Participant/(Feedback*technique)), data=pair_df)
emmeans <- emmeans(afmod, specs= "Feedback", by= "technique")
pairs(emmeans, adjust = "Tukey", infer = c(TRUE, TRUE))

emmeans1 <- emmeans(afmod, specs= "technique", by= "Feedback")
pairs(emmeans1, adjust = "Tukey", infer = c(TRUE, TRUE))

# pdf('Posthoc/Term_leak_diff_fb.pdf',  width = 6, height = 5)
# plot(emmeans, comparisons = TRUE, xlab=x)
# dev.off()




