# Author:                      Job van Riet
# Date:                        29-03-2022
# Function:                    Determine CABA-V7 patient characteristics.


# Import libraries ----

library(plyr)
library(dplyr)
library(ggplot2)
library(tableone)
library(patchwork)

# Helper theme.
source('R/misc_themes.R')

# Import data ----

data.Patient <- list()
data.Patient$Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$clinicalData <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical Characteristics') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$CAE <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Common Adverse Events') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$SAE <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Severe Adverse Events') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$Lab <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Lab') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$Hematology <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Hematology', col_types = c('text', 'date', 'numeric', 'numeric','numeric','numeric','numeric')) %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$PSA <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'PSA Measurements', col_types = c('text', 'date', 'numeric'))
data.Patient$priorChemo <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Prior Treatment - Chemo', col_types = c('text', 'text', 'date', 'date', 'numeric', 'text'))
data.Patient$priorHormonal <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Prior Treatment - Hormonal', col_types = c('text', 'text', 'date', 'date', 'text'))


# Table 1 - Patient Characteristics ----

# Add AR-V7 status (on baseline).
data.Patient$clinicalData <- data.Patient$clinicalData %>%
    # Add AR-V7 status (on baseline).
    dplyr::left_join(data.Patient$Overview %>% dplyr::distinct(`Subject Number`, `Inclusion (Treated with Caba)`, `AR-V7 (Baseline)`)) %>%
    dplyr::mutate(`AR-V7 (Baseline)` = ifelse(`Inclusion (Treated with Caba)` == 'Yes' & `AR-V7 (Baseline)` == 'Pos.', 'Main - Pos.', `AR-V7 (Baseline)`)) %>%

    dplyr::mutate(postTreatment_Clean = ifelse(grepl('Cabazitaxel|Enzalu|Abir', `Post-Treatment`), `Post-Treatment`, 'Other')) %>%
    dplyr::mutate(postTreatment_Clean = ifelse(`Post-Treatment` == 'Unknown', 'Unknown', postTreatment_Clean)) %>%
    dplyr::mutate(postTreatment_Clean = ifelse(grepl('Cabazitaxel', `Post-Treatment`), 'Cabazitaxel', postTreatment_Clean)) %>%

    # Determine which patients received prior local therapy
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(
        'Local therapy' = base::ifelse('Yes' %in% c(`Radical prostactomy`, `Bilateral orchiectomy`, `Prior internal RTX`, `Prior external radiotherapy prostate`), 'Yes', 'No')
    ) %>%
    dplyr::ungroup()

# Generate Table 1.
tableOne <- tableone::CreateTableOne(
    vars = c(
        'Age at registration',
        'WHO/ECOG PS at registration',
        'Total Gleason',
        'Local therapy',
        'Prior hormonal therapy',
        'Prior chemotherapy',
        'Prior RTX metastases',
        'Prior other treatment',
        'M-stage primary diagnosis',
        'PSA at primary diagnosis [ug/L]',
        'postTreatment_Clean'
    ),
    strata = c("AR-V7 (Baseline)"),
    data = data.Patient$clinicalData,
    smd = F,
    test = F,
    addOverall = T
)

# Print Table 1.
print(
    tableOne,
    test = F,
    nonnormal = c('Age at registration', 'Total Gleason', 'M-stage primary diagnosis', 'postTreatment_Clean', 'PSA at primary diagnosis [ug/L]'),
    exact = c('n'),
    smd = F,
    minMax = T,
    missing = T,
    explain = T,
    showAllLevels = T,
    quote = F,
    noSpaces = T
)


# Determine mean/median of lab-values ----

# PSA
data.Patient$clinicalData %>%
    dplyr::group_by(`AR-V7 (Baseline)`) %>%
    dplyr::summarise(
        mean = mean(na.omit(`PSA at primary diagnosis [ug/L]`)),
        sd = sd(na.omit(`PSA at primary diagnosis [ug/L]`)),
        out = sprintf('%s±%s', round(mean, 1), round(sd, 1))
    )

selectClosestDate <- function(data){
    data %>%
        dplyr::full_join(data.Patient$clinicalData %>% dplyr::select(`Subject Number`, `AR-V7 (Baseline)`, `Date: Registration`)) %>%
        dplyr::mutate(
            dateDiff = as.Date(as.character(`Date: Registration`, format="%Y/%m/%d")) - as.Date(as.character(`Date: Measurement`, format="%Y/%m/%d"))
        ) %>%
        dplyr::group_by(`Subject Number`) %>%
        dplyr::slice(
            ifelse(any(dateDiff >= 0), which.min(abs(dateDiff[dateDiff >= 0] - 0)), which.min(abs(dateDiff - 0)))
        ) %>%
        dplyr::ungroup() %>%
        dplyr::distinct() %>%
        dplyr::filter(dateDiff >= -60) %>%
        dplyr::full_join(data.Patient$clinicalData %>% dplyr::select(`Subject Number`, `AR-V7 (Baseline)`))
}

# Select measurement closest to date of registration.
dataHema <- selectClosestDate(data.Patient$Hematology %>% dplyr::filter(!is.na(HBG), !is.na(ANC)))
dataLab <- selectClosestDate(data.Patient$Lab %>% dplyr::filter(!is.na(ALP), !is.na(LDH)))

dataLab %>%
    dplyr::select(`Subject Number`, `AR-V7 (Baseline)`, ALP, LDH) %>%
    reshape2::melt() %>%
    dplyr::group_by(`AR-V7 (Baseline)`, variable) %>%
    dplyr::summarise(
        mean = mean(na.omit(value)),
        sd = sd(na.omit(value)),
        median = median(na.omit(value)),
        min = min(na.omit(value)),
        max = max(na.omit(value)),
        outMean = sprintf('%s±%s', round(mean, 1), round(sd, 1)),
        outMedian = sprintf('%s (%s - %s)', round(median, 1), round(min, 1), round(max, 1)),
        missing = sum(is.na(value))
    )


# Adverse events ----

data.Patient$CAE <- data.Patient$CAE %>% dplyr::filter(AEREL != 'Unrelated')
data.Patient$SAE <- data.Patient$SAE %>% dplyr::filter(SAEREL != 'Unrelated')

# Any adverse events.
dplyr::n_distinct(c(data.Patient$CAE %>% dplyr::pull(`Subject Number`), data.Patient$SAE %>% dplyr::pull(`Subject Number`)))

# Any adverse events >grade 3.
dplyr::n_distinct(c(data.Patient$CAE %>% dplyr::filter(AETOXGR >= 3) %>% dplyr::pull(`Subject Number`), data.Patient$SAE %>% dplyr::filter(SAETOXGR >= 3) %>% dplyr::pull(`Subject Number`)))

# Any serious events.
data.Patient$SAE %>% dplyr::summarise(totalN = dplyr::n_distinct(`Subject Number`))

# Events leading to treatment change.
data.Patient$SAE %>% dplyr::filter(SAEACN %in% c('Treatment interrupted', 'Treatment dose reduced')) %>% dplyr::summarise(totalN = dplyr::n_distinct(`Subject Number`))

# Common Adverse Events.
data.Patient$CAE %>% dplyr::group_by(AETERM) %>% dplyr::summarise(totalN = dplyr::n_distinct(`Subject Number`)) %>% dplyr::arrange(-totalN)

# Severe Adverse Events.
data.Patient$SAE %>% dplyr::group_by(SAETERM) %>% dplyr::summarise(totalN = dplyr::n_distinct(`Subject Number`)) %>% dplyr::arrange(-totalN)


# Determine CTC response ----

# CTC response is positive if CTC-count <5 in follow-up CTC count during caba-treatment.
responseCTC <- data.Patient$Overview %>%
    dplyr::filter(
        `Inclusion (Treated with Caba)` == 'Yes',
        (!is.na(`CTC – Treatment 3`) | !is.na(`CTC – Treatment 4`) | !is.na(`CTC – Treatment 5`))
    ) %>%
    dplyr::select(`Subject Number`, `CTC Count (Baseline – 7.5mL)`, `CTC – Treatment 3`, `CTC – Treatment 4`, `CTC – Treatment 5`) %>%
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(

        # Check if any CTC-measurement is below five CTC after caba-treatment.
        ctcResponse = ifelse(min(`CTC – Treatment 3`, `CTC – Treatment 4`, `CTC – Treatment 5`, na.rm = T) < 5, 'Formal Response (# CTC after any caba-treatment < 5)', 'No Response'),

        # Check if any CTC-measurement after caba-treatment is below 50% of the baseline CTC count.
        ctcDecline = ifelse(any((na.omit(c(`CTC – Treatment 3`, `CTC – Treatment 4`, `CTC – Treatment 5`)) / `CTC Count (Baseline – 7.5mL)`) <= .5), 'Formal Decline (50% drop of baseline CTC (7.5 mL) in any caba-treatment)', 'No Decline')

    ) %>%
    dplyr::ungroup()

# Figure 2 - Primary Endpoint ----

# Convert PSA measurements.
data.PSA <- data.Patient$PSA %>%

    # Filter on caba-treated patients.
    dplyr::full_join(data.Patient$clinicalData, by = 'Subject Number') %>%
    dplyr::full_join(data.Patient$Overview, by = 'Subject Number') %>%
    dplyr::filter(`Inclusion (Treated with Caba).x` == 'Yes') %>%

    # Convert dates.
    dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%

    # Center on start of caba. treatment.
    dplyr::mutate(
        dayPSA.Caba = `Date: PSA measurement` - `Date: Start of Post-Treatment`,
        dayDeath.Caba = `Date: Death` - `Date: Start of Post-Treatment`,

        # Correct for actual end of study.
        `Date: End of cabazitaxel cycles` = `Date: End of cabazitaxel cycles` + 21,

        duringCaba = ifelse(dayPSA.Caba > 0 & `Date: PSA measurement` <= `Date: End of cabazitaxel cycles`, 'Yes', 'No')
    ) %>%

    # Center on PSA prior to caba.
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(
        PSA.First = ifelse(any(dayPSA.Caba <= 0), tail(`PSA [ug/L]`[dayPSA.Caba <= 0], 1), -9999999999)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(`Subject Number`) %>%
    dplyr::mutate(
        PSA.Centered = (`PSA [ug/L]` / max(PSA.First)) - 1,
        dayEndOfCabaCycles = `Date: End of cabazitaxel cycles` - `Date: Start of Post-Treatment`,
        maxPSA.Decrease = `PSA [ug/L]`[duringCaba == 'Yes'][which.min(abs(`PSA [ug/L]`[duringCaba == 'Yes']))],
        maxDeltaPSA.Decrease = ifelse(all(duringCaba == 'No'), NA, maxPSA.Decrease - max(PSA.First))
    ) %>%
    dplyr::ungroup() %>%

    # Threshold PSA increase.
    dplyr::mutate(PSA.Centered = ifelse(PSA.Centered > 1, 1, PSA.Centered)) %>%

    # Sort patients on total length of caba-treatment.
    dplyr::mutate(`Subject Number` = factor(`Subject Number`, levels = unique(`Subject Number`[order(-dayEndOfCabaCycles)])))


# Remove data prior to CABA-V7 (except the first PSA) and post-CABA-V7
data.PSA <- data.PSA %>%
    dplyr::filter(duringCaba == 'Yes' | (`PSA [ug/L]` == PSA.First & dayPSA.Caba <= 0))

# PSA per caba. treatment (per patient, centered on start of caba.).
plot.PSA <- data.PSA %>%
    ggplot2::ggplot(aes(x = dayPSA.Caba, y = PSA.Centered)) +

    # Initial pre-caba PSA (baseline).
    ggplot2::geom_hline(yintercept = 0, color = '#F6A18C', alpha = .5, lty = 11, size = .25) +

    # Up/low line of CABA-V7 cycles.
    ggplot2::geom_rect(aes(ymin = -1.5, ymax = 1.5, xmin = 0, xmax = dayEndOfCabaCycles), fill = NA, color = 'royalblue2', lty = '11', lwd = .05) +

    # PSA measurements.
    ggplot2::geom_line(size = .4, lty = 11, color = 'grey50') +
    ggplot2::geom_point(size = .75, color = 'black') +

    # Show which PSA-point used as the baseline.
    ggplot2::geom_point(data = data.PSA %>% dplyr::filter(`PSA [ug/L]` == PSA.First, dayPSA.Caba <= 0), size = .9, color = 'royalblue2') +

    # Show which PSA-point used as the max. diff decrease PSA.
    ggplot2::geom_point(data = data.PSA %>% dplyr::filter(`PSA [ug/L]` == maxPSA.Decrease, duringCaba == 'Yes'), size = .9, color = '#E87E0D') +

    # Show days of death after Caba treatment.
    # ggplot2::geom_point(data = data.PSA %>% dplyr::filter(!is.na(`Date: Death`)) %>% dplyr::distinct(`Subject Number`, dayDeath.Caba), aes(x = dayDeath.Caba, y = 0), color = 'black', shape='\u2020', size = 3) +

    # Scales.
    ggplot2::scale_y_continuous(breaks = c(-1, 0, 1), labels = c('-100%', '0%', '≥100%'), expand = c(.1, .1)) +
    ggplot2::scale_color_discrete(guide = 'none') +
    scale_x_continuous(breaks = seq(-25, 301, 25), expand = c(0,0), limits = c(-30, 310)) +

    # Labs.
    ggplot2::labs(x = 'Number of days on cabazitaxel treatment', y = 'ΔPSA (in %), µg/L<br><span style = "font-size:6pt">(Centered on PSA prior to Cabazitaxel treatment)</span>') +

    # Facet + themes
    facet_grid(`Subject Number` ~ ., scales = 'free_x', switch = 'y') +
    theme_Job +
    theme(
        axis.text.y = element_text(size = 4.5),
        axis.text.x = element_text(size = 6, angle = 45),
        panel.spacing.y = ggplot2::unit(.1, "cm"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(linetype = 'dotted', color = '#D9D9D930'),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_rect(color = '#D9D9D9')
    )

# Show initial + delta PSA
plot.deltaPSA <- data.PSA %>%
    dplyr::select(`Subject Number`, 'Prior PSA (µg/L)' = PSA.First, 'Max. ΔPSA Decrease (%)' = maxDeltaPSA.Decrease) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
        `Max. ΔPSA Decrease (%)` = sprintf('(%s%s%%)', ifelse(`Max. ΔPSA Decrease (%)` >= 0, '↥', '↧'), round((`Max. ΔPSA Decrease (%)` / `Prior PSA (µg/L)`) * 100)),
        `Subject Number` = factor(`Subject Number`, levels = rev(levels(data.PSA$`Subject Number`)))
    ) %>%
    reshape2::melt(id.vars = 'Subject Number') %>%

    ggplot2::ggplot(aes(x = variable, y = `Subject Number`, label = value, fill = value,)) +

    # Add heatmap.
    ggplot2::geom_tile(width = 2, height = 1, colour = 'white', fill = 'white', lwd = .25, na.rm = T) +
    ggplot2::geom_text(size = 2.5, hjust = 0) +

    # Labs
    ggplot2::labs(x = NULL, y = NULL) +

    # Scales.
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +

    theme_Job +
    theme(
        text = element_text(family = 'Helvetica'),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank()
    )

# Show CTC / PSA response per patient.
plot.Response <- data.PSA %>%
    dplyr::distinct(`Subject Number`, `Response PSA`, `Response CTC`, `Response CTC-Decline`) %>%
    reshape2::melt(id.vars = 'Subject Number') %>%
    dplyr::mutate(
        `Subject Number` = factor(`Subject Number`, levels = rev(levels(data.PSA$`Subject Number`))),
        variable = factor(variable, levels = c('Response PSA', 'Response CTC', 'Response CTC-Decline'))
    ) %>%

    ggplot2::ggplot(aes(x = variable, y = `Subject Number`, fill = value)) +

    # Add heatmap.
    ggplot2::geom_tile(width = .75, height = .75, colour = 'grey25', lwd = .25, na.rm = T) +

    # Labs
    ggplot2::labs(x = NULL, y = NULL) +

    # Scales.
    ggplot2::scale_x_discrete(expand = c(0,0)) +
    ggplot2::scale_y_discrete(expand = c(0,0)) +
    ggplot2::scale_fill_manual(values = c('Formal Decline (50% drop of baseline CTC (7.5 mL) in any caba-treatment)' = '#659A37', 'Formal Response (# CTC after any caba-treatment < 5)' = '#659A37', 'Formal Response (Two consecutive PSA decreases ≤ 50% from baseline PSA)' = '#659A37', 'No Decline' = 'white', 'No Response' = 'white', '.' = 'grey75'), na.value = 'grey75', guide = guide_legend(title = 'Response Type', title.position = 'top', title.hjust = .5, nrow = 4)) +

    theme_Job +
    ggplot2::theme(
        text = element_text(family = 'Helvetica'),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank()
    )


# Combine plots -------------

## Combine landscape tracks.
layout <- "AAAB"

plot.PSA + plot.Response +
    patchwork::plot_layout(design = layout, heights = c(.4, .4, .4, 5, rep(.1, 8)), guides = 'collect') +
    patchwork::plot_annotation(tag_levels = 'a')


# Calculate CI --------------

# PSA response.
z <- binom.test(3, n = 25, 0.5, alternative="two.sided", conf.level = 0.95)
sprintf('PSA Response (3 of 25): p = %s (.95 CI: %s - %s)', z$p.value, round(z$conf.int[[1]], 3), round(z$conf.int[[2]], 3))

# CTC response.
z <- binom.test(3, n = 20, 0.5, alternative="two.sided", conf.level = 0.95)
sprintf('CTC Response (3 of 20): p = %s (.95 CI: %s - %s)', z$p.value, round(z$conf.int[[1]], 3), round(z$conf.int[[2]], 3))

# CTC decline.
z <- binom.test(10, n = 20, 0.5, alternative="two.sided", conf.level = 0.95)
sprintf('CTC Decline (10 of 20): p = %s (.95 CI: %s - %s)', z$p.value, round(z$conf.int[[1]], 3), round(z$conf.int[[2]], 3))


# Comparison to TROPIC ------

# Test the PSA-response rates.
test <- data.frame(
    row.names = c('CABA-V7', 'TRPOC'),
    PSA = c(3, 129),
    noPSA = c(22, 200)
)

chisq.test(test, simulate.p.value = T)
fisher.test(test, simulate.p.value = T, alternative = 'two.sided')