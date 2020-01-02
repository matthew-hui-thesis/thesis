list.of.packages <- c("shiny", "shinyjs", "dplyr", "data.table", "ggplot2", "tictoc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(shiny)
library(shinyjs)
library(dplyr)
library(data.table)
library(ggplot2)
library(tictoc)

#Input criteria
max.Ninf <- 12 #(minimum is 1)
max.Nobs <- 12 #(minimum is 2)

#Functions for EBE estimation
##Model parameters
TVCL=5.1659477035671495
TVV1=21.374461429274231
TVQ=3.15474007965699313E-002
TVV2=1.0023996233833952
B_HT_V1=2.3057725912876594
B_HT_CL=1.5870770112063479
B_WT_Q=0.74370418594815568
B_AGE_V2=0.65169699975111128
B_CLCR_CL=0.35717576675049739
B_DPA_CL=0.4072793470064626

SG=0.09787562
OMCL=3.60862468703888190E-003
OMV1=4.81269772393493056E-003
OMQ=3.46182573885426836E-002
OMV2=7.08896529053811425E-002
OMIOV=2.82142889480972402E-002

OM = diag(c(OMCL, OMV1, OMQ, OMV2))
SG = matrix(SG)

NsC1 = function(s, Rate, CL, V1, Q, V2, C10, C20) {
  V1*V2**2*C10*s**3 +
    (V1*V2*C10*Q+Rate*V2**2+V2**2*C20*Q+V1*V2*C10*Q)*s**2 +
    (2*Rate*V2*Q+V2*C20*Q**2+V1*C10*Q**2)*s +
    Rate*Q**2
}
comp2InfC1 <- function(t, Rate, popCL, etaCL, etaIOV, popV1, etaV1, popQ, etaQ, popV2, etaV2, C10, C20) {
  CL = popCL*exp(etaCL+etaIOV)
  V1 = popV1*exp(etaV1)
  Q = popQ*exp(etaQ)
  V2 = popV2*exp(etaV2)
  a = (Q/V2+Q/V1+CL/V1-sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  b = (Q/V2+Q/V1+CL/V1+sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  (Rate*Q/a/b-NsC1(-Q/V2, Rate, CL, V1, Q, V2, C10, C20)/(Q/V2)/(a-Q/V2)/(b-Q/V2)*exp(-Q/V2*t)+
      NsC1(-a, Rate, CL, V1, Q, V2, C10, C20)/a/(V2*a-Q)/(b-a)*exp(-a*t)+
      NsC1(-b, Rate, CL, V1, Q, V2, C10, C20)/b/(V2*b-Q)/(a-b)*exp(-b*t))/V1/V2
}
NsC2 = function(s, Rate, CL, V1, Q, V2, C10, C20) {
  V1*V2*C20*s**2 +
    (V2*Q*C20+V2*CL*C20+V1*Q*C10)*s +
    Rate*Q
}
comp2InfC2 <- function(t, Rate, popCL, etaCL, etaIOV, popV1, etaV1, popQ, etaQ, popV2, etaV2, C10, C20) {
  CL = popCL*exp(etaCL+etaIOV)
  V1 = popV1*exp(etaV1)
  Q = popQ*exp(etaQ)
  V2 = popV2*exp(etaV2)
  a = (Q/V2+Q/V1+CL/V1-sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  b = (Q/V2+Q/V1+CL/V1+sqrt((Q/V2+Q/V1+CL/V1)**2-4*CL*Q/V1/V2))/2
  (Rate*Q/a/b+NsC2(-a, Rate, CL, V1, Q, V2, C10, C20)/a/(a-b)*exp(-a*t)+
      NsC2(-b, Rate, CL, V1, Q, V2, C10, C20)/b/(b-a)*exp(-b*t))/V1/V2
}

getComp2conc <- function(times, timeInf, Rate, arrPopCL, arrPopV1, arrPopQ, arrPopV2, etaCL, etaV1, etaQ, etaV2) {
  #Next observation carried backward is assumed
  if (!any(times %in% timeInf)&&timeInf<max(times)) {
    added <- T
    iID <- (sort(c(times, timeInf))==timeInf) %>% which
    cTimes <- sort(c(times, timeInf))
    cRate <- c(rep(Rate, iID), rep(0, length(cTimes)-iID))
    cPopCL <- c(arrPopCL[1:(iID-1)], arrPopCL[iID], arrPopCL[iID:length(arrPopCL)])
    cPopV1 <- c(arrPopV1[1:(iID-1)], arrPopV1[iID], arrPopV1[iID:length(arrPopV1)])
    cPopQ <- c(arrPopQ[1:(iID-1)], arrPopQ[iID], arrPopQ[iID:length(arrPopQ)])
    cPopV2 <- c(arrPopV2[1:(iID-1)], arrPopV2[iID], arrPopV2[iID:length(arrPopV2)])
  } else {
    added <- F
    cTimes <- times
    cRate <- sapply(cTimes, function(t) {if (t<=timeInf) {Rate} else {0}})
    cPopCL <- arrPopCL
    cPopV1 <- arrPopV1
    cPopQ <- arrPopQ
    cPopV2 <- arrPopV2
  }
  arrC1 <- 0
  arrC2 <- 0
  for (i in 2:length(cTimes)) {
    lst_args <- list(t = cTimes[i]-cTimes[i-1], Rate = cRate[i],
                     popCL = cPopCL[i], etaCL = etaCL, etaIOV = 0,
                     popV1 = cPopV1[i], etaV1 = etaV1,
                     popQ = cPopQ[i], etaQ = etaQ,
                     popV2 = cPopV2[i], etaV2 = etaV2,
                     arrC1[i-1], arrC2[i-1])
    arrC1 <- append(arrC1, do.call("comp2InfC1", lst_args))
    arrC2 <- append(arrC2, do.call("comp2InfC2", lst_args))
  }
  if (added) {arrC1 <- arrC1[-iID]; arrC2 <- arrC2[-iID]}
  return(matrix(c(arrC1, arrC1), byrow = F, ncol = 2))
}

EBE_optim <- function(ETA, DATA) {
  with(DATA, {
    etaCL = ETA[1]
    etaV1 = ETA[2]
    etaQ = ETA[3]
    etaV2 = ETA[4]
    lapply(seq(Ninf), function(i) {
      lst_args = list(
        times = TAD[[i]],
        timeInf = infTime[i],
        Rate = dosePerBSA[i]*BSAdose[i]*2200.51/infTime[i],
        arrPopCL = TVCL*(HT[[i]]/156)**B_HT_CL*(CLCR[[i]]/179)**B_CLCR_CL*(dosePerBSA[i]/12)**B_DPA_CL,
        arrPopV1 = TVV1*(HT[[i]]/156)**B_HT_V1,
        arrPopQ = TVQ*(WT[[i]]/40)**B_WT_Q,
        arrPopV2 = TVV2*(AGE[[i]]/146)**B_AGE_V2,
        etaCL = etaCL, etaV1 = etaV1, etaQ = etaQ, etaV2 = etaV2
      )
      return(do.call(getComp2conc, lst_args)[-1,])
    }) %>% do.call(rbind, .) -> FGH
    F = FGH[,1] %>% matrix(., ncol = 1)
    H = FGH[,2] %>% matrix(., ncol = 1)
    COV <- ((H) %*% SG %*% t(H)) %>% diag %>% diag(., ncol=length(F), nrow=length(F))
    Y <- do.call(c, lapply(seq(Ninf), function(i) {DV[[i]][-1]}))
    RES <- Y-F
    OBJ <- log(det(COV)) + t(RES) %*% solve(COV) %*% RES + t(ETA) %*% solve(OM) %*% ETA
    return(OBJ %>% as.vector)
  })
}

#Function for dose optimization
AUC_optim <- function(DPA, DATA) {
  with(DATA, {
    P.LT.AUC_UB <- 1-((DPA*nowBSA*2200.51/ICL_MDPA/(DPA/12)**B_DPA_CL/AUC_UB) %>% log %>% pnorm(., 0, sqrt(OMIOV)))
    P.LT.AUC_LB <- 1-((DPA*nowBSA*2200.51/ICL_MDPA/(DPA/12)**B_DPA_CL/AUC_LB) %>% log %>% pnorm(., 0, sqrt(OMIOV)))
    P.IN <- P.LT.AUC_UB - P.LT.AUC_LB
    return(P.IN)
  })
}

MTX_UI <- fluidPage(
  shinyjs::useShinyjs(),
  #Hidden panel for example button use
    conditionalPanel("input.Ninf > 999.99", fluidRow(disabled(numericInput("egDummy", NULL, 0)))),
  #Hidden panel for result panel use
    conditionalPanel("input.Ninf > 999.99", fluidRow(disabled(numericInput("resDummy", NULL, 0)))),
  #Title
    fluidRow(column(6, titlePanel("Methotrexate dose calculator for Osteosarcoma"))),
    hr(),
  #Tab formatting
  tags$style(type = 'text/css', HTML(".container-fluid > .nav > li > a[data-value='tabOpt'] {background-color: #339966; color:white;}
                                       .container-fluid > .nav > li > a:hover[data-value='tabOpt'] {background-color: #228855; color:white;}
                                       .container-fluid > .nav > li > a:focus[data-value='tabOpt'] {background-color: #117744; color:white;}")),
  #tabs
  navbarPage(title = "", id = "selTab",
  tabPanel("Step 1 - Getting individual parameters", icon = icon("user", "fa-2x"), value = "tab1", wellPanel(
    #Getting individual parameters
    ##Number of infusions
    fluidRow(numericInput("Ninf", "Total number of infusions to input:", 1,
                          min = 1, max = max.Ninf, step = 1) %>% column(., width = 3),
             column(3, actionButton("exampleButton", "Example values (ID=11)"))),
    ##Use renderUI to produce the UI for the max.Ndose doses
    tags$div(id = "InfPanel")
  )),
  tabPanel("Step 2 - Setting therapeutic target(s)", icon = icon("bullseye", "fa-2x"), value = "tab2", wellPanel(
    ##Adjust therapeutic target
    fluidRow(column(12, textOutput("crteriaDetail"))),
    fluidRow(titlePanel(" ")),
    fluidRow(column(12, textOutput("tut_str_def_AUC"))),
    fluidRow(column(3, numericInput("def_AUC_LB", "Lower bound:", 4000, min = 0, step = 100)),
             column(3, numericInput("def_AUC_UB", "Upper bound:", 12000, min = 0, step = 100)))
  )),
  tabPanel("Step 3 - Current demographics", icon = icon("sliders", "fa-2x"), value = "tab3", wellPanel(
    #Info for this dose
    fluidRow(column(6,titlePanel("Current demographics"))),
    wellPanel(
      fluidRow(column(2, textInput("nowBSA", "BSA (m2):", 1)),
               column(2, textInput("nowHT", "HT (cm):", 1)),
               column(2, textInput("nowCLCR", "CLCR (mL/min/1.73m2):", 1)),
               column(2, textInput("nowWT", "WT (kg):", 1)),
               column(2, textInput("nowAGE", "AGE (months):", 1)))
    )
  )),
  tabPanel("DOSE OPTIMIZATION", icon = icon("cogs", "fa-2x"), value = "tabOpt", wellPanel(
    #Calculation button
    fluidRow(column(actionButton("calculate", "Calculate", style="color: #ffffff; background-color: #228855; border-color: #11b45d; font-size: 18px"), width = 2)),
    fluidRow(column(6, textOutput("msgCal"))),
    br(),
    br(),
    #Result part
    conditionalPanel("input.resDummy == 1",
                     wellPanel(fluidRow(column(6, titlePanel("Results"))),
                               fluidRow(column(12, wellPanel(style = {"background-color: #fff6ed; border-color: #514e4b"},
                                 fluidRow(column(6, fluidRow(column(12, textOutput("optTitle"))),
                                                 fluidRow(column(12, textOutput("optDPA"))),
                                                 fluidRow(column(12, textOutput("optDose"))),
                                                 hr(),
                                                 fluidRow(column(12, textOutput("probTitle"))),
                                                 fluidRow(column(12, textOutput("Dose_P_AUC_LB_95"))),
                                                 fluidRow(column(12, textOutput("Dose_P_AUC_UB_95"))),
                                                 hr(),
                                                 fluidRow(column(12, textOutput("selProbTitle"))),
                                                 fluidRow(column(12, textOutput("optDose_P_LT_AUC_LB"))),
                                                 fluidRow(column(12, textOutput("optDose_P_IN"))),
                                                 fluidRow(column(12, textOutput("optDose_P_GT_AUC_UB")))),
                                          column(6, fluidRow(plotOutput("plot_1")),
                                                 fluidRow(textOutput("Legend_1")))
                               )))),
                               hr(),
                               wellPanel(
                                 style = {"background-color: #eaf7ff; border-color: #3d382c"},
                                 fluidRow(column(12, textOutput("profTitle"))),
                                 fluidRow(column(12, plotOutput("plot_2")))
                               ),
                               hr()
                     )
    ),
    #Styling
    tags$head(tags$style("#optTitle{font-size: 22px; text-decoration: underline}")),
    tags$head(tags$style("#optDPA{font-size: 20px}")),
    tags$head(tags$style("#optDose{font-size: 20px}")),
    tags$head(tags$style("#probTitle{font-size: 22px; text-decoration: underline}")),
    tags$head(tags$style("#optDose_P_LT_AUC_LB{font-size: 18px}")),
    tags$head(tags$style("#optDose_P_IN{font-size: 18px}")),
    tags$head(tags$style("#optDose_P_GT_AUC_UB{font-size: 18px}")),
    tags$head(tags$style("#selProbTitle{font-size: 22px; text-decoration: underline}")),
    tags$head(tags$style("#Dose_P_AUC_LB_95{font-size: 18px}")),
    tags$head(tags$style("#Dose_P_AUC_UB_95{font-size: 18px}")),
    tags$head(tags$style("#profTitle{font-size: 22px; text-decoration: underline}")),
    tags$head(tags$style("#Legend_1{font-size: 16px}")),
    tags$head(tags$style("#crteriaDetail{font-size: 18px}")),
    tags$head(tags$style("#tut_str_def_AUC{font-size: 16px}"))
  )))
)

MTX_Server <- function(input, output, session) {
#Initiation
  withProgress(
    {incProgress(0, message = "Loading... just a few seconds...")
    },
    message = "Calculator setup"
  )
  output$crteriaDetail <- renderText(paste0("Dose is optimized by maximizing the probability of AUC within the therapeutic target of ", format(input$def_AUC_LB, nsmall = 0, big.mark = ","), " to ", format(input$def_AUC_UB, nsmall = 0, big.mark = ","), " uM*h."))
  output$tut_str_def_AUC <- renderText("Or you may define the AUC target yourself:")
#UI creation
  change_Ninf <- reactive({input[["Ninf"]]})
  ExistNinf <<- 0
  
  #Creating inf panels
    observeEvent(change_Ninf(), {
      Ninf <- change_Ninf()
      if (Ninf>ExistNinf) {    
        #Decide list to insert
          if (ExistNinf==0) {lstToIns <- seq(Ninf)} else {lstToIns <- seq(Ninf)[-seq(ExistNinf)]}
        #Update list inserted
          ExistNinf <<- Ninf
        #Insert UI & reactive component
          lapply(lstToIns,
                 function(InfN) {
            #Insert UI
              insertUI(
                selector = '#InfPanel',
                ui = conditionalPanel(paste0("input.Ninf >= ", InfN),
                                      hr(),
                                      wellPanel(
                                        style = if (InfN%%2==0) {"background-color: #ccf5fe; border-color: #44535c"} else {"background-color: #aed7f1; border-color: #44535c"},
                                        fluidRow(column(3, titlePanel(paste0("Dose #", InfN)))),
                                        fluidRow(column(3, numericInput(paste0("Nobs_", InfN),
                                                                        "Number of observations (including start of infusion):", 2, min = 2, max = max.Nobs, step = 1)),
                                                 column(3, textInput(paste0("dosePerBSA_", InfN), "Total dose per BSA (g/m2):", 12)),
                                                 column(3, textInput(paste0("BSA_", InfN), "BSA (m2):", 0.84)),
                                                 column(3, textInput(paste0("infTime_", InfN), "Infusion time (hr):", 6))
                                        ),
                                        fluidRow(column(2, "Time after infusion started (hr):"),
                                                 column(2, "Conc (umol/L):"),
                                                 column(1, "BSA (m2):"),
                                                 column(1, "HT (cm):"),
                                                 column(2, "CLCR (mL/min/1.73m2):"),
                                                 column(1, "WT (kg):"),
                                                 column(1, "AGE (months):")
                                        ),
                                        fluidRow(
                                          column(1, offset = 4, checkboxInput(paste0("fix_BSA_", InfN), "Same", FALSE)),
                                          column(1, checkboxInput(paste0("fix_HT_", InfN), "Same", FALSE)),
                                          column(2, checkboxInput(paste0("fix_CLCR_", InfN), "Same", FALSE)),
                                          column(1, checkboxInput(paste0("fix_WT_", InfN), "Same", FALSE)),
                                          column(1, checkboxInput(paste0("fix_AGE_", InfN), "Same", FALSE))
                                        ),
                                        tags$div(id = paste0("ObsPanel_", InfN))
                                      ))
              )
            #Update Obs reactive components
              assign(paste0("change_Nobs_", InfN), reactive({input[[paste0("Nobs_", InfN)]]}), inherits = T)
            #Creating Obs components
              if (!exists(paste0("ExistNobs_", InfN))) {assign(paste0("ExistNobs_", InfN), 0)}
              observeEvent(get(paste0("change_Nobs_", InfN))(), {
                Nobs <- get(paste0("change_Nobs_", InfN))()
                ExistNobs <- get(paste0("ExistNobs_", InfN))
                lstToIns_obs <- numeric()
                if (Nobs>ExistNobs) {
                  #Decide list to insert
                    if (ExistNobs==0) {lstToIns_obs <- seq(Nobs)} else {lstToIns_obs <- seq(Nobs)[-seq(ExistNobs)]}
                  #Update list inserted
                    assign(paste0("ExistNobs_", InfN), Nobs, inherits = T)
                  #Insert UI
                    lapply(lstToIns_obs, function(ObsN) {
                      if (ObsN==1) {
                        insertUI(
                          selector = paste0("#ObsPanel_", InfN),
                          ui = conditionalPanel(paste0("input.Nobs_", InfN, " >= ", 1),
                                                fluidRow(
                                                  column(2, disabled(textInput(paste0("TAD_", InfN, "_", 1), NULL, 0))),
                                                  column(2, disabled(textInput(paste0("DV_", InfN, "_", 1), NULL, 0))),
                                                  column(1, textInput(paste0("BSA_", InfN, "_", 1), NULL, 1)),
                                                  column(1, textInput(paste0("HT_", InfN, "_", 1), NULL, 1)),
                                                  column(2, textInput(paste0("CLCR_", InfN, "_", 1), NULL, 1)),
                                                  column(1, textInput(paste0("WT_", InfN, "_", 1), NULL, 1)),
                                                  column(1, textInput(paste0("AGE_", InfN, "_", 1), NULL, 1))
                                                )),
                          session = session
                        )
                      } else {
                        insertUI(
                          selector = paste0("#ObsPanel_", InfN),
                          ui = conditionalPanel(paste0("input.Nobs_", InfN, " >= ", ObsN),
                                                fluidRow(
                                                  column(2, textInput(paste0("TAD_", InfN, "_", ObsN), NULL,
                                                                      (if (ObsN==2) {6} else {(ObsN-2)*24}))),
                                                  column(2, textInput(paste0("DV_", InfN, "_", ObsN), NULL, 1)),
                                                  column(1, {if (input[[paste0("fix_BSA_", InfN)]]) {
                                                    disabled(textInput(paste0("BSA_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("BSA_", InfN, "_", ObsN), NULL, 1)
                                                  }}),
                                                  column(1, {if (input[[paste0("fix_HT_", InfN)]]) {
                                                    disabled(textInput(paste0("HT_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("HT_", InfN, "_", ObsN), NULL, 1)
                                                  }}),
                                                  column(2, {if (input[[paste0("fix_CLCR_", InfN)]]) {
                                                    disabled(textInput(paste0("CLCR_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("CLCR_", InfN, "_", ObsN), NULL, 1)
                                                  }}),
                                                  column(1, {if (input[[paste0("fix_WT_", InfN)]]) {
                                                    disabled(textInput(paste0("WT_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("WT_", InfN, "_", ObsN), NULL, 1)
                                                  }}),
                                                  column(1, {if (input[[paste0("fix_AGE_", InfN)]]) {
                                                    disabled(textInput(paste0("AGE_", InfN, "_", ObsN), NULL, 1))
                                                  } else {
                                                    textInput(paste0("AGE_", InfN, "_", ObsN), NULL, 1)
                                                  }})
                                                ))
                        )
                      }
                    })
                    lapply(
                      lstToIns_obs,
                      function(ObsN) {
                        if (ObsN!=1) {
                          assign(paste0("react_fix_BSA_", InfN, "_", ObsN), reactive({input[[paste0("fix_BSA_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_BSA_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("BSA_", InfN, "_", ObsN), !(get(paste0("react_fix_BSA_", InfN, "_", ObsN))()))
                          })
                          assign(paste0("react_fix_HT_", InfN, "_", ObsN), reactive({input[[paste0("fix_HT_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_HT_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("HT_", InfN, "_", ObsN), !(get(paste0("react_fix_HT_", InfN, "_", ObsN))()))
                          })
                          assign(paste0("react_fix_CLCR_", InfN, "_", ObsN), reactive({input[[paste0("fix_CLCR_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_CLCR_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("CLCR_", InfN, "_", ObsN), !(get(paste0("react_fix_CLCR_", InfN, "_", ObsN))()))
                          })
                          assign(paste0("react_fix_WT_", InfN, "_", ObsN), reactive({input[[paste0("fix_WT_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_WT_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("WT_", InfN, "_", ObsN), !(get(paste0("react_fix_WT_", InfN, "_", ObsN))()))
                          })
                          assign(paste0("react_fix_AGE_", InfN, "_", ObsN), reactive({input[[paste0("fix_AGE_", InfN)]]}), inherits = T)
                          observeEvent(get(paste0("react_fix_AGE_", InfN, "_", ObsN))(), {
                            shinyjs::toggleState(paste0("AGE_", InfN, "_", ObsN), !(get(paste0("react_fix_AGE_", InfN, "_", ObsN))()))
                          })
                        }
                      }
                    )
                }
              })
          })
      }
  })
    
    
    
#Example button
  observeEvent(input$exampleButton, {
    withProgress({
      updateTextInput(session, "Ninf", value = 2)
      updateTextInput(session, "nowBSA", value = 2.05)
      updateTextInput(session, "nowHT", value = 185)
      updateTextInput(session, "nowCLCR", value = 168.35)
      updateTextInput(session, "nowWT", value = 83.3)
      updateTextInput(session, "nowAGE", value = 175)
      updateNumericInput(session, "egDummy", value = 1)
    }, message = "Inputinge example values")
  })
  ch_egDummy <- reactive({input$egDummy})
  observeEvent(ch_egDummy(), {
    if (ch_egDummy()==1) {
      withProgress({
        updateTextInput(session, "Nobs_1", value = 5)
        updateNumericInput(session, "egDummy", value = 2)
      }, message = "Inputinge example values")
    } else if (ch_egDummy()==2) {
      withProgress({
        updateTextInput(session, "Nobs_2", value = 5)
        updateNumericInput(session, "egDummy", value = 3)
      }, message = "Inputinge example values")
    } else if (ch_egDummy()==3) {
      withProgress({
        updateCheckboxInput(session, "fix_BSA_1", value = TRUE)
        updateCheckboxInput(session, "fix_HT_1", value = TRUE)
        updateCheckboxInput(session, "fix_CLCR_1", value = FALSE)
        updateCheckboxInput(session, "fix_WT_1", value = TRUE)
        updateCheckboxInput(session, "fix_AGE_1", value = TRUE)
        updateCheckboxInput(session, "fix_BSA_2", value = TRUE)
        updateCheckboxInput(session, "fix_HT_2", value = TRUE)
        updateCheckboxInput(session, "fix_CLCR_2", value = FALSE)
        updateCheckboxInput(session, "fix_WT_2", value = TRUE)
        updateCheckboxInput(session, "fix_AGE_2", value = TRUE)
        updateNumericInput(session, "egDummy", value = 4)
      }, message = "Inputinge example values")
    } else if (ch_egDummy()==4) {  
      withProgress({
        updateTextInput(session, "dosePerBSA_1", value = 12)
        updateTextInput(session, "BSA_1", value = 2.05)
        updateTextInput(session, "infTime_1", value = 6)
        
        updateTextInput(session, "TAD_1_2", value = 6)
        updateTextInput(session, "TAD_1_3", value = 24)
        updateTextInput(session, "TAD_1_4", value = 48)
        updateTextInput(session, "TAD_1_5", value = 72)
        updateTextInput(session, "DV_1_2", value = 600)
        updateTextInput(session, "DV_1_3", value = 3.7)
        updateTextInput(session, "DV_1_4", value = 0.38)
        updateTextInput(session, "DV_1_5", value = 0.11)
        updateTextInput(session, "BSA_1_1", value = 2.05)
        updateTextInput(session, "HT_1_1", value = 185)
        updateTextInput(session, "CLCR_1_1", value = 158.9972)
        updateTextInput(session, "CLCR_1_2", value = 158.9972)
        updateTextInput(session, "CLCR_1_3", value = 158.9972)
        updateTextInput(session, "CLCR_1_4", value = 194.0305)
        updateTextInput(session, "CLCR_1_5", value = 170.8627)
        updateTextInput(session, "WT_1_1", value = 83.3)
        updateTextInput(session, "AGE_1_1", value = 175)
        
        updateTextInput(session, "dosePerBSA_2", value = 13)
        updateTextInput(session, "BSA_2", value = 2.05)
        updateTextInput(session, "infTime_2", value = 6)
        
        updateTextInput(session, "TAD_2_2", value = 6)
        updateTextInput(session, "TAD_2_3", value = 24)
        updateTextInput(session, "TAD_2_4", value = 48)
        updateTextInput(session, "TAD_2_5", value = 72)
        updateTextInput(session, "DV_2_2", value = 730)
        updateTextInput(session, "DV_2_3", value = 10)
        updateTextInput(session, "DV_2_4", value = 0.79)
        updateTextInput(session, "DV_2_5", value = 0.21)
        updateTextInput(session, "BSA_2_1", value = 2.05)
        updateTextInput(session, "HT_2_1", value = 185)
        updateTextInput(session, "CLCR_2_1", value = 168.35)
        updateTextInput(session, "CLCR_2_2", value = 168.35)
        updateTextInput(session, "CLCR_2_3", value = 124.4326)
        updateTextInput(session, "CLCR_2_4", value = 150.6289)
        updateTextInput(session, "CLCR_2_5", value = 154.7)
        updateTextInput(session, "WT_2_1", value = 83.3)
        updateTextInput(session, "AGE_2_1", value = 175)
        updateNumericInput(session, "egDummy", value = 0)
      }, message = "Inputinge example values")
    }
  })

#Calculate button
  observeEvent(input$calculate, {
    updateNumericInput(session, "resDummy", value = 0)
    ErrFree <- TRUE
    withProgress({
      tic("Done reading values!")
      tryCatch({
  #Reading input values
          incProgress(1/4, message = "Reading input values")
          as.object <- function(name) {
            eval(parse(text = paste0("input$", name)))
          }
        #AUC
          AUC_LB <- as.object("def_AUC_LB") %>% as.numeric
          AUC_UB <- as.object("def_AUC_UB") %>% as.numeric
          if (is.numeric(AUC_LB)&&(AUC_LB>0)) {} else {stop()}
          if (is.numeric(AUC_UB)&&(AUC_UB>0)&&(AUC_UB>AUC_LB)) {} else {stop()}
        #Ninf
          if ((input$Ninf%%1==0)&&(input$Ninf>0)) {Ninf <- input$Ninf} else {stop()}
        #Nobs, dosePerBSA, BSA, infTime
          Nobs <- integer(Ninf)
          dosePerBSA <- numeric(Ninf)
          BSAdose <- numeric(Ninf)
          infTime <- numeric(Ninf)
          TAD <- list()
          DV <- list()
          BSA <- list()
          HT <- list()
          CLCR <- list()
          WT <- list()
          AGE <- list()
          for (i in 1:Ninf) {
            Nobs[i] <- as.object(paste0("Nobs_", i))
            dosePerBSA[i] <- as.object(paste0("dosePerBSA_", i)) %>% as.numeric
            BSAdose[i] <- as.object(paste0("BSA_", i)) %>% as.numeric
            infTime[i] <- as.object(paste0("infTime_", i)) %>% as.numeric
            if ((Nobs[i]%%1==0)&&(Nobs[i]>1)) {} else {stop()}
            if (is.numeric(dosePerBSA[i])&&(dosePerBSA[i]>0)) {} else {stop()}
            if (is.numeric(BSAdose[i])&&(BSAdose[i]>0)) {} else {stop()}
            if (is.numeric(infTime[i])&&(infTime[i]>0)) {} else {stop()}
        #TAD, DV, BSA, HT, CLCR, WT, AGE
            curTAD <- numeric(Nobs[i])
            curDV <- numeric(Nobs[i])
            curBSA <- numeric(Nobs[i])
            curHT <- numeric(Nobs[i])
            curCLCR <- numeric(Nobs[i])
            curWT <- numeric(Nobs[i])
            curAGE <- numeric(Nobs[i])
            for (j in 1:Nobs[i]) {
              curTAD[j] <- as.object(paste0("TAD_", i, "_", j)) %>% as.numeric
              curDV[j] <- as.object(paste0("DV_", i, "_", j)) %>% as.numeric
              curBSA[j] <- as.object(paste0("BSA_", i, "_", j)) %>% as.numeric
              curHT[j] <- as.object(paste0("HT_", i, "_", j)) %>% as.numeric
              curCLCR[j] <- as.object(paste0("CLCR_", i, "_", j)) %>% as.numeric
              curWT[j] <- as.object(paste0("WT_", i, "_", j)) %>% as.numeric
              curAGE[j] <- as.object(paste0("AGE_", i, "_", j)) %>% as.numeric
              if (as.object(paste0("fix_BSA_", i))) {curBSA[j] <- curBSA[1]}
              if (as.object(paste0("fix_HT_", i))) {curHT[j] <- curHT[1]}
              if (as.object(paste0("fix_CLCR_", i))) {curCLCR[j] <- curCLCR[1]}
              if (as.object(paste0("fix_WT_", i))) {curWT[j] <- curWT[1]}
              if (as.object(paste0("fix_AGE_", i))) {curAGE[j] <- curAGE[1]}
              if (is.numeric(curTAD[j])&&(((j==1)&&(curTAD[j]==0))||((j!=1)&&(curTAD[j]>0)))) {} else {stop()}
              if (is.numeric(curDV[j])&&(((j==1)&&(curDV[j]==0))||((j!=1)&&(curDV[j]>0)))) {} else {stop()}
              if (is.numeric(curBSA[j])&&(curBSA[j]>0)) {} else {stop()}
              if (is.numeric(curHT[j])&&(curHT[j]>0)) {} else {stop()}
              if (is.numeric(curCLCR[j])&&(curCLCR[j]>0)) {} else {stop()}
              if (is.numeric(curWT[j])&&(curWT[j]>0)) {} else {stop()}
              if (is.numeric(curAGE[j])&&(curAGE[j]>0)) {} else {stop()}
            }
            if (identical(curTAD, sort(curTAD))&&((curTAD %>% unique %>% length)==(curTAD %>% length))) {} else {stop()}
            TAD <- c(TAD, list(curTAD))
            DV <- c(DV, list(curDV))
            BSA <- c(BSA, list(curBSA))
            HT <- c(HT, list(curHT))
            CLCR <- c(CLCR, list(curCLCR))
            WT <- c(WT, list(curWT))
            AGE <- c(AGE, list(curAGE))
          }
        #BSA, HT, CLCR, WT, AGE (for this dose)
          nowBSA <- as.object("nowBSA") %>% as.numeric
          nowHT <- as.object("nowHT") %>% as.numeric
          nowCLCR <- as.object("nowCLCR") %>% as.numeric
          nowWT <- as.object("nowWT") %>% as.numeric
          nowAGE <- as.object("nowAGE") %>% as.numeric
          if (is.numeric(nowBSA)&&(nowBSA>0)) {} else {stop()}
          if (is.numeric(nowHT)&&(nowHT>0)) {} else {stop()}
          if (is.numeric(nowCLCR)&&(nowCLCR>0)) {} else {stop()}
          if (is.numeric(nowWT)&&(nowWT>0)) {} else {stop()}
          if (is.numeric(nowAGE)&&(nowAGE>0)) {} else {stop()}
        #ending message
          output$msgCal <- renderText("Inputs values okay!")
    toc()
  #Getting individual parameters
    tic("Obtaining individual estimates")
    incProgress(1/4, message = "Estimating individual parameters")
    resOptim <- optim(par = c(0,0,0,0), fn = EBE_optim, method = "L-BFGS-B",
                      DATA = list(
                        Ninf = Ninf, TAD = TAD, infTime = infTime,
                        dosePerBSA = dosePerBSA, BSAdose = BSAdose,
                        HT = HT, BSA = BSA, CLCR = CLCR, WT = WT, AGE = AGE,
                        DV = DV
                      ))
    assign("inspectOptim", resOptim, envir = .GlobalEnv)
    #Get etas
      ETA_CL <- resOptim$par[1]
      ETA_V1 <- resOptim$par[2]
      ETA_Q <- resOptim$par[3]
      ETA_V2 <- resOptim$par[4]
      ICL_MDPA <- TVCL*(nowHT/156)**B_HT_CL*(nowCLCR/179)**B_CLCR_CL*exp(ETA_CL)
      IV1 <- TVV1*(nowHT/156)**B_HT_V1*exp(ETA_V1)
      IQ <- TVQ*(nowWT/40)**B_WT_Q*exp(ETA_Q)
      IV2 <- TVV2*(nowAGE/146)**B_AGE_V2*exp(ETA_V2)
    toc()
  #Optimizing dose
    tic("Optimizing dose")
    incProgress(1/4, message = "Optimizing dose")
      #Optimal dose
        resOptim_DPA <- optim(par = 12, fn = AUC_optim, gr = NULL, method = "L-BFGS-B", lower = 1, upper = 50,
                              control = list(fnscale = -1),
                              DATA = list(nowBSA = nowBSA, ICL_MDPA = ICL_MDPA, AUC_UB = AUC_UB, AUC_LB = AUC_LB))
        final_optDose <- round(resOptim_DPA$par[1], 1)
        final_P.LT.AUC_LB <-
          round(1-((final_optDose*nowBSA*2200.51/ICL_MDPA/(final_optDose/12)**B_DPA_CL/AUC_LB) %>%
                     log %>% pnorm(., 0, sqrt(OMIOV))), 4)
        final_P.GT.AUC_UB <-
          round((final_optDose*nowBSA*2200.51/ICL_MDPA/(final_optDose/12)**B_DPA_CL/AUC_UB) %>%
                  log %>% pnorm(., 0, sqrt(OMIOV)), 4)
        final_P.IN <- 1 - final_P.GT.AUC_UB - final_P.LT.AUC_LB
      #Getting P95 dose information
        dose_P95_AUC_LB <- (ICL_MDPA*AUC_LB*(1/12)**B_DPA_CL*exp(qnorm(0.975,0,sqrt(OMIOV)))/2200.51/nowBSA)**(1/(1-B_DPA_CL))
        dose_P95_AUC_UB <- (ICL_MDPA*AUC_UB*(1/12)**B_DPA_CL*exp(qnorm(0.025,0,sqrt(OMIOV)))/2200.51/nowBSA)**(1/(1-B_DPA_CL))
      #Getting expected PK profile at optimal dose (by simulation)
        timePt <- seq(0.1, 96, by = 0.1)
        #Predicted expected Cp
          ICL <- ICL_MDPA*(final_optDose/12)**B_DPA_CL
          Ik12 <- IQ/IV1
          Ik21 <- IQ/IV2
          Ik10 <- ICL/IV1
          Ibeta <- (Ik12+Ik21+Ik10-sqrt((Ik12+Ik21+Ik10)**2-4*Ik21*Ik10))/2
          Ialpha <- Ik21*Ik10/Ibeta
          IA <- (Ialpha-Ik21)/IV1/(Ialpha-Ibeta)
          IB <- (Ibeta-Ik21)/IV1/(Ibeta-Ialpha)
          ExpectedConc <- sapply(timePt, function(t) {
            if (t<=6) {
              (final_optDose*2200.51*nowBSA/6*(IA/Ialpha*(1-exp(-Ialpha*t))+IB/Ibeta*(1-exp(-Ibeta*t)))) %>% return()
            } else {
              (final_optDose*2200.51*nowBSA/6*(IA/Ialpha*(1-exp(-Ialpha*6))*exp(-Ialpha*(t-6))+IB/Ibeta*(1-exp(-Ibeta*6))*exp(-Ibeta*(t-6)))) %>% return()
            }
          })
          PredProf <- data.frame(timePt, ExpectedConc)
          sampleSize <- 5000
          ETA_IOV <- rnorm(sampleSize, 0, sqrt(OMIOV))
          EPS <- rnorm(sampleSize*10, 0, sqrt(SG))
          PredConc <- sapply(ETA_IOV, function(eta_IOV) {
            ICL <- ICL_MDPA*(final_optDose/12)**B_DPA_CL*exp(eta_IOV)
            Ik12 <- IQ/IV1
            Ik21 <- IQ/IV2
            Ik10 <- ICL/IV1
            Ibeta <- (Ik12+Ik21+Ik10-sqrt((Ik12+Ik21+Ik10)**2-4*Ik21*Ik10))/2
            Ialpha <- Ik21*Ik10/Ibeta
            IA <- (Ialpha-Ik21)/IV1/(Ialpha-Ibeta)
            IB <- (Ibeta-Ik21)/IV1/(Ibeta-Ialpha)
            CpFcn <- function(t) {
              (final_optDose*2200.51*nowBSA/6*(IA/Ialpha*(1-exp(-Ialpha*6))*exp(-Ialpha*(t-6))+IB/Ibeta*(1-exp(-Ibeta*6))*exp(-Ibeta*(t-6)))) %>% return()
            }
            Conc_24h <- CpFcn(24)*exp(sample(EPS, 1))
            Conc_48h <- CpFcn(48)*exp(sample(EPS, 1))
            Conc_72h <- CpFcn(72)*exp(sample(EPS, 1))
            return(c(Conc_24h, Conc_48h, Conc_72h))
          }) %>% t %>% as.array %>% c(rep(24,sampleSize), rep(48,sampleSize), rep(72, sampleSize), .) %>% matrix(., ncol = 2) %>% as.data.frame
          colnames(PredConc) <- c("timePt", "Conc")
      #Generating graphs
        #AUC-dose graph
          if (ceiling(dose_P95_AUC_LB - 5)==0) {dose_LB <- 1} else {dose_LB <- max(ceiling(dose_P95_AUC_LB - 5), 0.1)}
          dose_UB <- ceiling(dose_P95_AUC_UB + 3)
          dose_seq <- seq(dose_LB, dose_UB, 0.1)
          AUC_P025 <- dose_seq*2200.51*nowBSA/ICL_MDPA/(dose_seq/12)**B_DPA_CL/exp(qnorm(0.975,0,sqrt(OMIOV)))
          AUC_Mid <- dose_seq*2200.51*nowBSA/ICL_MDPA/(dose_seq/12)**B_DPA_CL
          AUC_P975 <- dose_seq*2200.51*nowBSA/ICL_MDPA/(dose_seq/12)**B_DPA_CL/exp(qnorm(0.025,0,sqrt(OMIOV)))
          plot_1 <- ggplot(data.frame(dose_seq, AUC_P025, AUC_Mid, AUC_P975), aes(x=dose_seq, y=AUC_Mid)) +
                      geom_line(aes(y=AUC_Mid), linetype = "solid", size = 1) +
                      geom_hline(yintercept=AUC_LB, linetype = "dashed", size = 1) +
                      geom_segment(x=round(dose_P95_AUC_LB, 1), xend=round(dose_P95_AUC_LB, 1),
                                   y=0, yend=AUC_LB, linetype = "dotted", size = 1) +
                      geom_segment(x=final_optDose, xend=final_optDose,
                                   y=0, yend=AUC_Mid[which(dose_seq==final_optDose)], linetype = "dotted", size = 1) +
                      geom_segment(x=round(dose_P95_AUC_UB, 1), xend=round(dose_P95_AUC_UB, 1),
                                   y=0, yend=AUC_UB, linetype = "dotted", size = 1) +
                      geom_hline(yintercept=AUC_UB, linetype = "dashed", size = 1) +
                      geom_ribbon(aes(ymin=AUC_P025, ymax=AUC_P975), fill="blue", alpha=0.25) +
                      geom_ribbon(aes(ymin=AUC_LB, ymax=AUC_UB), fill="green", alpha=0.05) +
                      scale_x_continuous(name = "Dose (g/m2)",
                                         limits = c(min(dose_seq), max(dose_seq)),
                                         breaks = c(min(dose_seq), round(dose_P95_AUC_LB, 1), final_optDose,
                                                    round(dose_P95_AUC_UB, 1), max(dose_seq))) +
                      scale_y_continuous(name = "AUC (uM*h)",
                                         limits = c(min(AUC_P025), max(AUC_P975)),
                                         breaks = c(AUC_LB, (AUC_LB+AUC_UB)/2, AUC_UB)) +
                      ggtitle("Predicted AUC against dose per BSA") +
                      theme(axis.text=element_text(size=18),
                            axis.title=element_text(size=21,face="bold"),
                            title = element_text(size = 16, vjust = 0.5),
                            panel.background = element_rect(fill = "#fff6ed",colour = NA),
                            panel.grid.minor = element_blank(),
                            panel.grid.major = element_blank(),
                            plot.background = element_rect(fill = "#fff6ed",colour = NA))

          
          plot_2 <- ggplot(PredProf, aes(x=timePt, y=ExpectedConc)) +
            stat_boxplot(data = PredConc, aes(group = factor(timePt), y = Conc), geom = "errorbar", width = 10) +
            geom_boxplot(data = PredConc, aes(group = factor(timePt), y = Conc), width = 10) + 
            geom_line(size = 1) + 
            scale_x_continuous(name = "Time after dose (hours)", breaks = c(6,24,48,72,96)) +
            scale_y_continuous(name = "Predicted plasma concentration (uM/L)", trans = "log10", breaks = c(0.1,1,10,100)) + 
            ggtitle("Predicted plasma concentration profile") +
            theme(axis.text=element_text(size=18),
                  axis.title.x = element_text(size=21,face="bold"),
                  axis.title.y = element_text(size=18,face="bold"),
                  title = element_text(size=16, vjust = 0.5),
                  panel.background = element_rect(fill = "#eaf7ff",colour = NA),
                  plot.background = element_rect(fill = "#eaf7ff",colour = NA))
      toc()
  #Generating output
    tic("Generating output")
    incProgress(1/4, message = "Generating output")
      output$optTitle <- renderText("Optimal dose")
      output$optDPA <- renderText(paste0("Optimal dose per BSA = ", final_optDose, " g/m2"))
      output$optDose <- renderText(paste0("Optimal dose (rounded to 2 d.p.) = ", round((final_optDose * nowBSA),2), " g"))
      output$probTitle <- renderText("Dose range with highest probability")
      output$optDose_P_LT_AUC_LB <- renderText(paste0("Probability of less than ", format(AUC_LB, nsmall = 0, big.mark = ","), " uM*h = ", final_P.LT.AUC_LB * 100, "%"))
      output$optDose_P_IN <- renderText(paste0("Probability of between ", format(AUC_LB, nsmall = 0, big.mark = ","), "~", format(AUC_UB, nsmall = 0, big.mark = ","), " uM*h = ", final_P.IN * 100, "%"))
      output$optDose_P_GT_AUC_UB <- renderText(paste0("Probability of greater than ", format(AUC_UB, nsmall = 0, big.mark = ","), " uM*h = ", final_P.GT.AUC_UB * 100, "%"))
      output$selProbTitle <- renderText("Probabilities at optimal dose")
      output$Dose_P_AUC_LB_95 <- renderText(paste0("Minimum dose to ensure that P(AUC>", format(AUC_LB, nsmall = 0, big.mark = ","), ") > 95% = ",
                                                 round(dose_P95_AUC_LB, 1), " g/m2",
                                                 " (or ", round((dose_P95_AUC_LB * nowBSA),2), " g)"))
      output$Dose_P_AUC_UB_95 <- renderText(paste0("Maximum dose to ensure that P(AUC<", format(AUC_UB, nsmall = 0, big.mark = ","), ") > 95% = ",
                                                  round(dose_P95_AUC_UB, 1), " g/m2",
                                                 " (or ", round((dose_P95_AUC_UB * nowBSA),2), " g)"))
      output$plot_1 <- renderPlot({plot_1})
      output$profTitle <- renderText("Expected profiles at optimal dose")
      output$plot_2 <- renderPlot({plot_2})
      output$Legend_1 <- renderText("Blue region shows the 95% confidence interval. Green region shows the desired therapeutic target.")
    toc()
    output$msgCal <- renderText("Done!")
      }, warning = function(w) {
      }, error = function(e) {
        output$msgCal <- renderText("Problems with input values! Please check!")
      }, finally = {
      })
    })
    updateNumericInput(session, "resDummy", value = 1)
  })
}

shinyApp(MTX_UI, MTX_Server)