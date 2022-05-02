# Author:                      Job van Riet
# Date:                        29-03-2022
# Function:                    Generate Survival Curves for the CABA-V7 study.

# Import libraries ----

library(plyr)
library(dplyr)
library(ggplot2)
library(survminer)
library(gtsummary)
library(extrafont)
library(patchwork)


# Functions ----

# Helper theme.
source('R/misc_themes.R')

# Generate survival plots with p-values and median OS.
plotSurvival <- function(fit, ylim, data, palette = 'jco', hr = NULL){
    
    # Generate survival plot.
    x <- survminer::ggsurvplot(
        fit = fit,
        pval = F,
        size = .825,
        break.time.by = 10,
        break.y.by = .2,
        palette = palette,
        risk.table = T,
        tables.height = .3,
        xlab = 'Time (in months)',
        axes.offset = F,
        ylim = c(0, 1.05),
        xlim = c(0, ylim), strata.labels = 'c',
        risk.table.col = 'strata', censor.shape = '+',
        fontsize = 3,
        risk.table.title = 'No. at risk',
        ggtheme = ggplot2::theme(
            legend.position = 'bottom',
            legend.direction = 'horizontal',
            text = ggplot2::element_text(size=9, family='Helvetica', face = 'bold'),
            axis.title.x = ggtext::element_textbox_simple(width = NULL, halign = .5),
            axis.title.y = ggtext::element_textbox_simple(size = 8, orientation = 'left-rotated', width = NULL, halign = .5),
            panel.grid.major.x = ggplot2::element_line(colour = 'grey90', linetype = 'dotted'),
            panel.grid.major.y = ggplot2::element_line(colour = '#E5E5E5', linetype = 'dotted'),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.background = ggplot2::element_rect(fill = NA, colour = 'black'),
            legend.text = ggtext::element_markdown()
        )
    )
    
    # Add the log-rank p-value.
    p.logrank <- survminer::surv_pvalue(fit = fit, method = 'log-rank', data = data, test.for.trend = F)
    x$plot <- x$plot + ggplot2::annotate("text", x = max(x$data.survplot$time) * .75, y = 1, label = paste0('log-rank: ', p.logrank$pval.txt), size = 2.5)
    
    # Add HR (if two groups)
    if(!is.null(hr)){
        
        HR.CI <- round(summary(hr)$conf.int, 2)
        HR.p <- round(summary(hr)$waldtest[[3]], 2)
        HR.CI <- sprintf('HR (.95%% CI): %s (%s - %s)', HR.CI[[1]], HR.CI[[3]], HR.CI[[4]])
        x$plot <- x$plot + ggplot2::annotate("text", x = max(x$data.survplot$time) * .75, y = .9, label = HR.CI, size = 2.5)
    }
    
    # Add the median OS.
    medianOS <- x$data.survplot %>%
        dplyr::group_by(strata) %>%
        dplyr::summarise(
            medianOS = round(median(time, na.rm = T), 2),
            label = sprintf('%s - %s mo.', unique(strata), medianOS)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(-medianOS)
    
    x$plot <- x$plot + ggplot2::annotate("text", x = max(x$data.survplot$time) * .75, y = .75, label = paste0('Median OS (Desc.):\n', paste(medianOS$label, collapse = '\n')), size = 2.5)
    
    # Remove legends.
    x$plot <- x$plot + theme_Job + theme(legend.position = 'none')
    x$table <- x$table + theme_Job + theme(legend.position = 'none')
    
    return(x)
}

plotHR <- function(data, withQ = F){
    x <- data %>% 
        gtsummary::tbl_regression(
            exponentiate = T, 
            add_estimate_to_reference_rows = T,
        ) %>%
        gtsummary::add_n() %>% 
        gtsummary::add_global_p() %>% 
        gtsummary::add_nevent() %>% 
        gtsummary::bold_p() %>% 
        gtsummary::bold_labels() %>% 
        gtsummary::italicize_levels() %>% 
        gtsummary::sort_p() %>% 
        bstfun::add_inline_forest_plot(header = '', spec_pointrange.args = list(lim = c(-3, 3), width = 550, cex = 1, col = "black", pch = 1))
    
    if(withQ){
        x <- x %>% gtsummary::add_q()
    }
    
    x %>% bstfun::as_ggplot()
    
}

# Import data ----

data.Patient <- list()
data.Patient$Overview <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, skip = 1, sheet = 'Sample Overview') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))
data.Patient$clinicalData <- readxl::read_excel('Misc./Suppl. Table 1 - Overview of Data.xlsx', trim_ws = T, sheet = 'Clinical Characteristics') %>% dplyr::mutate(`Subject Number` = as.character(`Subject Number`))


# Convert and clean data ----

data.Survival <- data.Patient$clinicalData %>%
    # Convert dates.
    dplyr::mutate_at(dplyr::vars(dplyr::contains('Date:')), as.Date) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        dateCensor = ifelse(!is.na(`Date: Death`), `Date: Death`, na.omit(c(`Date: Last follow-up`, `Date: End of study`, `Date: Pre-screening`))),
        dateCensor = as.Date(dateCensor, origin = '1970-01-01'),
        daysFromPreScreeningToEnd = dateCensor - `Date: Pre-screening`,
        monthsFromPreScreeningToEnd = daysFromPreScreeningToEnd / (365.25 / 12)
    ) %>%
    dplyr::inner_join(data.Patient$Overview) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
        `Dichotomized CTC count (Baseline)` = ifelse(`CTC Count (Baseline – 7.5mL)` >= 5, 'CTC Count (Baseline) ≥5', 'CTC Count (Baseline) <5'),
        `Patient recieved Cabazitaxel` = ifelse(grepl('Cabazitaxel', `Post-Treatment`), 'Yes', 'No'),
        `Inclusion (Treated with Cabazitaxel)` = ifelse(is.na(`Inclusion (Treated with Caba)`), 'No', `Inclusion (Treated with Caba)`),
        `WHO status (Pooled)` = ifelse(`WHO/ECOG PS at registration` %in% c(1,2), '1 - 2', `WHO/ECOG PS at registration`)
    )

plotFits <- list()


# Multivariate Cox-regression ---------------------------------------------

## Determine relevant factors (p <= 0.1; n = 125 complete cases) ----

data.AIC.All <- data.frame(data.Survival) %>%
    dplyr::select(
        monthsFromPreScreeningToEnd,
        Survival,
        AR.V7..Baseline.,
        WHO.ECOG.PS.at.registration,
        Max..VAF,
        Total.Gleason,
        Age.at.registration,
        Nr..of.coding.mutations,
        Dichotomized.CTC.count..Baseline.,
        WHO.status..Pooled.,
        Patient.recieved.Cabazitaxel
    ) %>%
    dplyr::filter(!is.na(Survival), !is.na(Max..VAF), !is.na(Total.Gleason))

MASS::stepAIC(survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline. + Patient.recieved.Cabazitaxel + Patient.recieved.Cabazitaxel*AR.V7..Baseline. + Max..VAF + Total.Gleason + Age.at.registration + Nr..of.coding.mutations + Dichotomized.CTC.count..Baseline. + WHO.status..Pooled., data = data.AIC.All, ties = 'breslow'), direction = 'backward')

# Pos. / Neg. Caba.
data.AIC.Caba <- data.frame(data.Survival) %>%
    dplyr::filter(
        Patient.recieved.Cabazitaxel == 'Yes',
        AR.V7..Baseline. != 'Und.'
    ) %>% 
    dplyr::select(
        monthsFromPreScreeningToEnd,
        Survival,
        AR.V7..Baseline.,
        WHO.ECOG.PS.at.registration,
        Max..VAF,
        Total.Gleason,
        Age.at.registration,
        Nr..of.coding.mutations,
        Dichotomized.CTC.count..Baseline.,
        WHO.status..Pooled.
    ) %>%
    dplyr::filter(!is.na(Survival), !is.na(Max..VAF), !is.na(Total.Gleason))

MASS::stepAIC(survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline. + Max..VAF + Total.Gleason + Age.at.registration + Nr..of.coding.mutations + Dichotomized.CTC.count..Baseline. + WHO.status..Pooled., data = data.AIC.Caba, ties = 'breslow') , direction = 'backward')


# Multivariate Cox Regression on relevant factors -------------------------

plotFits$multiCox.AllSamples <- data.Survival %>% 
    dplyr::filter(!is.na(Survival)) %>% 
    survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `AR-V7 (Baseline)` + `Patient recieved Cabazitaxel` + `AR-V7 (Baseline)` * `Patient recieved Cabazitaxel` + `Dichotomized CTC count (Baseline)`, data = ., ties = 'breslow') %>% 
    plotHR(., withQ = T)

plotFits$multiCox.PosNegCaba <- data.Survival %>%
    dplyr::filter(!is.na(Survival)) %>% 
    dplyr::filter(`Patient recieved Cabazitaxel` == 'Yes', `AR-V7 (Baseline)` != 'Und.') %>% 
    survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ `Dichotomized CTC count (Baseline)` + `AR-V7 (Baseline)`, data = ., ties = 'breslow') %>% 
    plotHR(., withQ = T)


# Survival Analysis (Cox regression ----

## Survival - AR-V7 (All included; n = 137) ----

fit.AllInClusion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.frame(data.Survival))
names(fit.AllInClusion$strata) <-  c('AR-V7 Neg.', 'AR-V7 Pos.', 'AR-V7 Und.')
plotFits$AllInclusion <- plotSurvival(fit.AllInClusion, data = data.frame(data.Survival), ylim = 51, palette = c('#648FFF', '#FE6100', '#4D4D4D'))

fit.AllInClusion.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count..Baseline., data = data.frame(data.Survival))
names(fit.AllInClusion.CTC$strata) <-  c('CTC Count (Baseline) <5', 'CTC Count (Baseline) ≥5')
plotFits$fit.AllInClusion.CTC <- plotSurvival(fit.AllInClusion.CTC, data = data.frame(data.Survival), ylim = 51, palette = c('#f23005', '#ffbe73'))


## Survival - AR-V7 (All included; n = 68) ----

data.CabaOnly <- data.Survival %>% 
    dplyr::filter(
        `AR-V7 (Baseline)` != 'Und.',
        `Patient recieved Cabazitaxel` == 'Yes'
    ) %>% 
    data.frame()

fit.CabaOnly <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline., data = data.CabaOnly)
names(fit.CabaOnly$strata) <-  c('AR-V7 Neg.', 'AR-V7 Pos.')
plotFits$CabaOnly <- plotSurvival(fit.CabaOnly, data = data.CabaOnly, ylim = 51, palette = c('#648FFF', '#FE6100'))

fit.CabaOnly.CTC <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Dichotomized.CTC.count..Baseline., data = data.CabaOnly)
names(fit.CabaOnly.CTC$strata) <-   c('CTC Count (Baseline) <5', 'CTC Count (Baseline) ≥5')
plotFits$fit.CabaOnly.CTC <- plotSurvival(fit.CabaOnly.CTC, data = data.CabaOnly, ylim = 51, palette = c('#f23005', '#ffbe73'))


# Generate Multi-Kaplan Figure (AR-V7) - All inclusion ----

layout <- "AA
BC
DE"

## All inclusion ----
plotFits$multiCox.AllSamples +
    plotFits$AllInclusion$plot +
    plotFits$fit.AllInClusion.CTC$plot +
    plotFits$AllInclusion$table +
    plotFits$fit.AllInClusion.CTC$table +
    patchwork::plot_layout(design = layout, heights = c(1.5, 1, .2), guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Helvetica'))

## Caba only Pos./Neg. ----
plotFits$multiCox.PosNegCaba +
    plotFits$CabaOnly$plot +
    plotFits$fit.CabaOnly.CTC$plot +
    plotFits$CabaOnly$table +
    plotFits$fit.CabaOnly.CTC$table +
    patchwork::plot_layout(design = layout, heights = c(1.5, 1, .2), guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Helvetica'))



# Suppl. OS ---------------------------------------------------------------

## WHO - All samples ----
fit.AllInClusion.WHO <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = data.frame(data.Survival))
names(fit.AllInClusion.WHO$strata) <-  c('WHO status: 0', 'WHO status: 1-2')
fit.AllInClusion.WHO.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ WHO.status..Pooled., data = data.frame(data.Survival))
plotFits$fit.AllInClusion.WHO <- plotSurvival(fit.AllInClusion.WHO, data = data.frame(data.Survival), hr = fit.AllInClusion.WHO.hr, ylim = 51, palette = c('#2C3D4F', '#1ABB9A'))

## CABA treatment  - All samples ----
fit.AllInClusion.Treatment <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Patient.recieved.Cabazitaxel, data = data.frame(data.Survival))
names(fit.AllInClusion.Treatment$strata) <-  c('No cabazitaxel treatment', 'With cabazitaxel treatment')
fit.AllInClusion.Treatment.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ Patient.recieved.Cabazitaxel, data = data.frame(data.Survival))
plotFits$fit.AllInClusion.Treatment <- plotSurvival(fit.AllInClusion.Treatment, data = data.frame(data.Survival), hr = fit.AllInClusion.Treatment.hr, ylim = 51, palette = c('blue', 'red'))

## Caba treatment vs. AR-V7 ----
fit.CabavsAR <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7..Baseline. + Patient.recieved.Cabazitaxel, data = data.Survival %>%  dplyr::filter(`AR-V7 (Baseline)` != 'Und.') %>% data.frame())
names(fit.CabavsAR$strata) <-  c('AR-V7 Neg. - No cabazitaxel treatment', 'AR-V7 Neg. - With cabazitaxel treatment', 'AR-V7 Pos. - No cabazitaxel treatment', 'AR-V7 Pos. - With cabazitaxel treatment')
plotFits$fit.CabavsAR <- plotSurvival(fit.CabavsAR, data = data.Survival %>%  dplyr::filter(`AR-V7 (Baseline)` != 'Und.') %>% data.frame(), ylim = 51, palette = c('#648FFF', 'darkblue', 'orange', '#FE6100'))

## Converters. ----
data.BetweenConversion <- data.Survival %>% dplyr::filter(`Inclusion (Treated with Cabazitaxel)` == 'Yes')
fit.BetweenConversion <- survminer::surv_fit(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.frame(data.BetweenConversion))
names(fit.BetweenConversion$strata) <-  c('AR-V7 Conversion (Pos. -> Neg.)', 'AR-V7 Retainment (Pos. -> Pos.)')
fit.BetweenConversion.hr <- survival::coxph(formula = survival::Surv(monthsFromPreScreeningToEnd, Survival) ~ AR.V7.Conversion, data = data.frame(data.BetweenConversion))
plotFits$fit.BetweenConversion <- plotSurvival(fit.BetweenConversion, data = data.BetweenConversion %>% data.frame(), ylim = 21, hr = fit.BetweenConversion.hr, palette = c('#00A94D', '#FE6100'))

### Combine Kaplans ----

layout <- "AB
EF
CD
GH"

plotFits$fit.AllInClusion.WHO$plot +
    plotFits$fit.AllInClusion.Treatment$plot +
    plotFits$fit.CabavsAR$plot +
    plotFits$fit.BetweenConversion$plot +
    plotFits$fit.AllInClusion.WHO$table +
    plotFits$fit.AllInClusion.Treatment$table +
    plotFits$fit.CabavsAR$table +
    plotFits$fit.BetweenConversion$table +
    patchwork::plot_layout(design = layout, heights = c(1, .2), guides = 'auto') +
    patchwork::plot_annotation(tag_levels = 'a') & ggplot2::theme(plot.tag = element_text(size = 11, family = 'Helvetica'))
