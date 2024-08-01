library(tidyverse)
library(ggpubr)
library(rstatix)
library(dplyr)
library(car)
library(rlang)

set.seed(123)


method2_wiltest <- function(df, x){
  condition_var= x
  counts =  df %>% group_by(Feedback) %>% 
    count(Feedback) %>% 
    summarise(min(n))
  
  # get data the same length
  pre_df = df %>%
    group_by(Feedback) %>% 
    select(condition_var = {{x}}, Feedback) %>%
    sample_n(min(counts$`min(n)`))
  
  pre_df = as.data.frame(pre_df)
  # Compute t-test
  pwc <- pre_df  %>% 
    # select(condition_var = {{x}}, technique, Feedback)%>%
    wilcox_test(condition_var ~Feedback, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance()  
  
  # Pairwise t-test
  pwc1 <- pre_df %>%
    pairwise_t_test(
      condition_var ~ Feedback, paired = TRUE, detailed = TRUE,
      p.adjust.method = "bonferroni"
    )
  
  # Summary results
  summary_stats <- pre_df%>%
    group_by(Feedback) %>%
    get_summary_stats(condition_var, type = "full")
  
  #Shapiro-Wilk normality test
  d <- with(pre_df, 
            condition_var[Feedback == "NF"] - condition_var[Feedback == "F"])
  # pre_df_sw <- pre_df %>%
  #   group_by(Feedback) %>%
  #   sample_n(5000)
  tryCatch(
    sw <- shapiro.test(d),

    error = function(ex){
      warning("Had more than 5000 entries - ", ex)
      #Shapiro-Wilk normality test for the differences
      sw <- shapiro.test(with(pre_df_sw, 
                              condition_var[Feedback == "NF"] - condition_var[Feedback == "F"]))
      })
      sw_df <- data.frame(statistic = sw$method,p_value = sw$p.value)
  
  
  # Visualisation
  pre_df %>%
    group_by(Feedback) %>%
    get_summary_stats(condition_var, type = "mean_sd")
  
  bxp <- ggboxplot(
    pre_df, x = "Feedback", y = "condition_var",
    add = "point",
    color = "Feedback",
    # main = paste(condition, x),
    palette = "jco"
  )
  bxp+theme(legend.position = "right",
            plot.margin = margin(1, 5, 1, 1))
  
  
  # Visualization: box plots with p-values
  # Prepare file path
  
  file_name <- sprintf("Plot/Pre-term/%s_wilcoxon.pdf", x)
  pwc <- pwc %>% add_xy_position(x = "Feedback")
  pdf(file_name,  width = 8, height = 5)
  print(
    bxp +
      stat_pvalue_manual(pwc,  label = "{p.adj} {p.adj.signif}") +
      ylab(condition_var) +
      labs(
        caption = get_pwc_label(pwc)
      )
  )
  dev.off()
  
  file_name <- sprintf("Plot/Pre-term/%s_ttest.pdf", x)
  pwc1 <- pwc1 %>% add_xy_position(x = "Feedback")
  pdf(file_name,  width = 8, height = 5)
  print(
    bxp +
      stat_pvalue_manual(pwc1,  label = "{p.adj} {p.adj.signif}") +
      ylab(condition_var) +
      labs(
        caption = get_pwc_label(pwc1)
      )
  )
  dev.off()
  return (list(summary_stats, pwc, pwc1, sw_df))
}


# Data preparation
df <- read.csv("Data/processed_data_final.csv")
# Select term or pre-term or all mannikin
mannikin =1
df <- df[df$mannikin == mannikin, ]

# "Ventilatory_rate"
cols_to_summarize = c("leak_diff_per", "Inspiratory_time", 
                      "Expiratory_time","Vti")
# Define empty data frame to store output

output_df_stat <- data.frame()
output_df_wilcoxon <- data.frame()
output_df_pwc <- data.frame()


# res = method2_wiltest(df,x)
# All case
for (x in cols_to_summarize){
    res = method2_wiltest(df, x)
    summary_stats = cbind(x, as.data.frame(res[1]))
    pwc = cbind(x, as.data.frame(res[2])[c("estimate", "group1","group2", 
                                           "n1","n2","p.adj")], res[4])
    pwc1 = cbind(x, as.data.frame(res[3])[c("estimate", "group1","group2", 
                                            "n1","n2","p.adj")], res[4])
    
    output_df_stat <- rbind(output_df_stat, summary_stats)
    output_df_wilcoxon <- rbind(output_df_wilcoxon, pwc)
    output_df_pwc <- rbind(output_df_pwc, pwc1)
    # output_df_sw <- rbind(output_df_sw, sw)
}


write.csv(output_df_stat, file = "Result/Pre-term_summary_stat.csv", row.names = FALSE)
write.csv(output_df_pwc, file = "Result/Pre-term_ttest.csv", row.names = FALSE)
write.csv(output_df_wilcoxon, file = "Result/Pre-term_wilcoxontest.csv", row.names = FALSE)