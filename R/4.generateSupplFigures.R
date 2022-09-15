# Author:                      Job van Riet
# Date:                        29-03-2022
# Function:                    Generate supplementary figures on the CABA-V7 cohort.


# Import libraries --------------------------------------------------------

library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)

# Helper theme.
source('R/misc_themes.R')


# Import data -------------------------------------------------------------

overviewPatients <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 1)

# Add number of patients per group.
overviewPatients <- overviewPatients %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>%
    dplyr::mutate(
        `AR-V7 (Baseline)` = ifelse(`AR-V7 (Baseline)` == 'Pos.', 'AR-V7<sup>Pos.</sup>', ifelse(`AR-V7 (Baseline)` == 'Neg.', 'AR-V7<sup>Neg.</sup>', 'AR-V7<sup>Und.</sup>')),
        `AR-V7 (Baseline) with n` = sprintf('%s<br>(<i>n</i> = %s)', `AR-V7 (Baseline)`, dplyr::n_distinct(`Subject Number`))
    ) %>%
    dplyr::ungroup()


# Correlation AR-V7 vs. CTC ----

# Do statistical tests.
stat.test <- overviewPatients %>%
    dplyr::mutate(g = `AR-V7 (Baseline) with n`, value = `CTC Count (Baseline)`) %>%
    rstatix::pairwise_wilcox_test(value ~ g, exact = T, p.adjust.method = 'none', detailed = T, paired = F) %>%
    rstatix::adjust_pvalue(method = 'BH') %>%
    rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

overviewPatients %>%
    dplyr::group_by(`AR-V7 (Baseline) with n`) %>%
    dplyr::mutate(medianCTC = median(`CTC Count (Baseline)`)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(., aes(x = reorder(`AR-V7 (Baseline) with n`, -medianCTC), y = `CTC Count (Baseline)`, fill = `AR-V7 (Baseline) with n`, label = medianCTC, group = `AR-V7 (Baseline) with n`)) +
    
    # Add half-half plots with median labels.
    gghalves::geom_half_boxplot(side = 'l', outlier.shape = NA, notch = F, show.legend = F) +
    gghalves::geom_half_point_panel(side = 'r', shape = 21, size = 1, color = 'black', position = ggbeeswarm::position_quasirandom(width = .15)) +
    ggplot2::stat_summary(fun=median, colour='black', geom='text', size = 3, show.legend = FALSE, vjust=-4, angle = 90, hjust = .5) +
    
    ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = stat.test, y.position = 9, step.increase = .02, tip.length = .01) +
    ggplot2::scale_fill_manual(values = c('#648FFF', '#FE6100', '#4D4D4D'), guide = 'none') +
    ggplot2::labs(x = 'AR-V7 determination (Baseline)', y = 'CTC Count (Baseline)') +
    ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), limits = c(-.25, 15000), breaks = c(0, 3, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000), expand = c(0,0)) +
    theme_Job
