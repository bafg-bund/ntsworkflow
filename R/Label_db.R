

# functions Label_DB ####

# Includes functions to create json-files for storaging Data and insert Data to the Label_DB



#' create_DB_json
#'
#' @param sampleList from NT-App/ environment
#' @param json_path path to json-files on Z
#' @param db_path path to both DBs on Server (Spectra_DB = db_v7 and Label_DB = db_v5)
#' @param grouped alignmentTable from NT-App/ environment
#' @param peaklist from NT-App/ environment
#' @param annotationTable from NT-App/ environment
#'
#' @return save Json for Label_DB from NT-App Data
#' @export
#' @import DBI
#' @import RSQLite
#'
#' @examples
create_DB_json <- function(samplelist_i=sampleList,
                          json_path="~/Data/Toni/Lable_Datenbank/Json_Files",
                          db_path="~/sqlite_local",
                          grouped_i=grouped,
                          peaklist_i=peaklist,
                          annotationTable_i=annotationTable,
                          datenlist_i=datenList){
  
  cat(
   "Die Importierung von Daten in die Label_DB erfolgt über Dateien im Format JSON
    Die Grundlage sind prozessierte Daten aus Kevins NTS-App, d.h. Schritte wie Blankabzug
    oder weitergehende Filter sollten an dieser Stelle bereits erledigt sein.
     \n")
  
  
  ## Aufbau der Grundstruktur ##
  # compGroupComp ist die many-to-many Verbindung zwischen compound und compoundGroup 
  compound <- data.frame(compound_id=integer(0), 
                         CAS=character(0), 
                         formula=character(0), 
                         SMILES=character(0), 
                         name=character(0), 
                         chem_list_id=character(0),
                         stringsAsFactors=FALSE) 
  
  compoundGroup <- data.frame(compoundGroup_id=c(1:16), 
                              name=as.character(c("BfG", "Pharmaceutical", "Transformation_product", "Antimicrobial", "Food_additive", "Fungicide",
                                                  "Herbicide", "Industrial_process", "Insecticide", "Metabolite", "Natural_product",
                                                  "Personal_care_product", "Pesticide", "LfU", "Pigment", "Surrogate_standard")),
                              stringsAsFactors=FALSE) 
  
  compGroupComp <- data.frame(compoundGroup_id=integer(0),
                              compound_id=integer(0))
  
  experiment <- data.frame(experiment_id=integer(0), 
                           compound_id=integer(0), 
                           parameter_id=integer(0),
                           sample_id=integer(0),
                           mz=double(0),
                           rt=double(0),
                           time_added=character(0), 
                           adduct=character(0), 
                           isotope=character(0),
                           Spectra_DB_score=integer(0),
                           stringsAsFactors=FALSE) 
  
  experimentGroup <- data.frame(experimentGroup_id=c(1,2), 
                                name=as.character(c("BfG", "LfU")),
                                stringsAsFactors=FALSE) 
  
  # expGroupExp ist die many-to-many Verbindung zwischen experiment und experimentGroup  
  expGroupExp <- data.frame(experimentGroup_id=integer(0), 
                            experiment_id=integer(0))
  
  fragment <- data.frame(fragment_id=integer(0), 
                         mz=double(0), 
                         int=double(0), 
                         experiment_id=integer(0)) 
  
  parameter <- data.frame(parameter_id=c(1:4), 
                          chrom_method=as.character(rep("dx.doi.org/10.1016/j.chroma.2015.11.014", times=4)),
                          instrument=as.character(c("LC-ESI-QTOF TripleTOF 5600 SCIEX", "LC-ESI-QTOF TripleTOF 5600 SCIEX",
                                                    "LC-ESI-QTOF TripleTOF 6600 SCIEX", "LC-ESI-QTOF TripleTOF 6600 SCIEX")), 
                          polarity=as.character(c("pos","neg","pos","neg")), 
                          ionisation=as.character(c("ESI","ESI","ESI","ESI")), 
                          CE=c(40,40,40,40), 
                          CES=c(15,15,15,15), 
                          ce_unit=as.character(c("V","V","V","V")), 
                          col_type=as.character(c("Q","Q","Q","Q")),
                          stringsAsFactors=FALSE) 
  
  sample <- data.frame(sample_id=integer(0), 
                       samplename=character(0), 
                       sample_type=character(0), 
                       location=character(0), 
                       enrichment=logical(0), 
                       enrichment_factor=integer(0), 
                       data_name=character(0), 
                       data_location=character(0), 
                       date=character(0), 
                       contact=character(0), 
                       project=character(0), 
                       leaching_id=integer(0), 
                       extraction_id=integer(0),
                       stringsAsFactors=FALSE)
  
  leaching <- data.frame(leaching_id=c(1), 
                         leaching_method=as.character(c("method XY")), 
                         sample_mass=NA, 
                         sample_area=NA,
                         stringsAsFactors=FALSE) 
  
  # sampLeachSamp ist die many-to-many Verbindung zwischen sample und leaching 
  sampLeachSamp <- data.frame(sample_id=integer(0), 
                              leaching_id=integer(0))
  
  extraction <- data.frame(extraction_id=c(1), 
                           extraction_method=as.character(c("dx.doi.org/10.1016/j.watres.2019.115366")),
                           sample_mass=NA,
                           stringsAsFactors=FALSE) 
  
  # samplExtractSamp ist die many-to-many Verbindung zwischen sample und extraction
  samplExtractSamp <- data.frame(sample_id=integer(0), 
                                 extraction_id=integer(0))
  
  ## Abfragen von Informationen über die Probe/ Messung ##  
  
  cat("\nWelche Proben sollen in einer JSON vereinigt werden?\n")
  show_samplelist <- data.frame("Number"=c(1:nrow(samplelist_i)),
                                "Filename"=stringr::str_match(basename(samplelist_i$File), "(.*)\\.mzXML")[,2])
  print(show_samplelist)
  cat("\nWelche Proben sollen in einer JSON vereinigt werden?")
  sample_i_start <- as.numeric(readline("Bitte Nummer der ersten Proben angeben ->   "))
  while (!any(sample_i_start==c(1:nrow(samplelist_i))) | is.na(sample_i_start)){
    cat("Nummer der ersten Probe ist außerhalb der Indizierung!")
    sample_i_start <- as.numeric(readline("Bitte Nummer der ersten Proben angeben ->   "))
  }
  sample_i_end <- as.numeric(readline("Bitte Nummer der letzten Proben angeben ->   "))
  while (!any(sample_i_end==c(1:nrow(samplelist_i))) | is.na(sample_i_end)){
    cat("Nummer der letzten Probe ist außerhalb der Indizierung!")
    sample_i_end <- as.numeric(readline("Bitte Nummer der letzten Proben angeben ->   "))
  }
  samples_i <- c(sample_i_start:sample_i_end)
  
  #  cat("\nIst die Probe in Replikaten erstellt/ gemessen worden?
  #    [1] nein      (n=1)
  #    [2] Duplikat  (n=2)
  #    [3] Triplikat (n=3)
  #    PS: Es wird davon ausgegangen, dass Replikate in einem Block prozessiert wurden!\n")
  #  repikat_i <- as.numeric(readline("Bitte zutreffende ID in die Console eingeben und bestätigen ->   "))
  #  if(!any(repikat_i==(c(1:3))) || repikat_i==c("") || repikat_i > length(peaklist_i) ){
  #    stop("Angegebene Anzahl von Replikaten ist nicht Teil der Einstellungsmöglichkeiten.\n
  #       Bitte Wert überprüfen oder ggf Skript anpassen!")}
  
  cat("\nWie sind die Messparameter?
    [1] Chromatographie: Gudrun, MS: LC-ESI-QTOF TripleTOF 5600 SCIEX - pos - ESI - CE=40 V - CES=15 V - Col. type= Q
    [2] Chromatographie: Gudrun, MS: LC-ESI-QTOF TripleTOF 5600 SCIEX - neg - ESI - CE=40 V - CES=15 V - Col. type= Q
    [3] Chromatographie: Gudrun, MS: LC-ESI-QTOF TripleTOF 6600 SCIEX - pos - ESI - CE=40 V - CES=15 V - Col. type= Q
    [4] Chromatographie: Gudrun, MS: LC-ESI-QTOF TripleTOF 6600 SCIEX - neg - ESI - CE=40 V - CES=15 V - Col. type= Q\n")
  parameter_id_i <- as.numeric(readline("Bitte zutreffende ID in die Console eingeben und bestätigen ->   "))
  while (!any(parameter_id_i==(c(1:4))) | is.na(parameter_id_i)){
    cat("\nDie angegebene ID steht nicht zur Auswahl!\n")
    parameter_id_i <- as.numeric(readline("Bitte zu zutreffende ID ZWISCHEN 1 UND 4 WÄHLEN! ->   "))
  }
    
  cat("\nWoher stammen die Proben?
    [1] BfG
    [2] LfU\n")
  experimentGroup_id_i <- as.numeric(readline("Bitte zutreffende ID in die Console eingeben und bestätigen ->   "))
  while (!any(experimentGroup_id_i==(c(1:2))) | is.na(experimentGroup_id_i)){
    cat("\nDie angegebene ID steht nicht zur Auswahl!\n")
    experimentGroup_id_i <- as.numeric(readline("Bitte zu zutreffende ID ZWISCHEN 1 UND 2 WÄHLEN! ->   "))
  }

  cat("\nUnter welchem Namen soll die Probe in der DB gespeichert werden?\n")
  samplename_i <- readline("Bitte Probenname in die Console eingeben und bestätigen ->   ")
  while (file.exists(paste0(json_path,"/",samplename_i,".json")) | samplename_i=="") {
    if (file.exists(paste0(json_path,"/",samplename_i,".json"))) {cat("\nEine Datei dieses Namens existiert bereits!\n")}
    if (samplename_i=="") {cat("\nFür die Speicherung der json braucht die Probe einen Namen!\n")}
    samplename_i <- readline("Bitte einen anderen Probenname angeben ->   ")
  }
  
  cat("\nIst die Probe aus einer Extraktion oder einem Leaching entstanden?
    [1] Extraktion
    [2] Leaching\n")
  sample_type_ii <- as.numeric(readline("Bitte zutreffende ID in die Console eingeben und bestätigen ->   "))
  while (!any(sample_type_ii==(c(1:2))) | is.na(sample_type_ii)){
    cat("\nDie angegebene ID steht nicht zur Auswahl!\n")
    sample_type_ii <- as.numeric(readline("Bitte zu zutreffende ID ZWISCHEN 1 UND 2 WÄHLEN! ->   "))
  }
  
  if (sample_type_ii==1) {
    sample_type_i <- "Extraction"
    cat("\nWelche Extraktionsmethode wurde verwendet?
    [1] Solvent1= MeOH/MQ - Solvent2= MeOH/2% FA - T=100 °C - p=100 bar (Lise Methode)\n")
    extraction_id_i <- as.numeric(readline("Bitte zutreffende ID in die Console eingeben und bestätigen ->   "))
    while (!any(extraction_id_i==(c(1:1))) | is.na(extraction_id_i)){
      cat("\nDie angegebene ID steht nicht zur Auswahl!\n")
      extraction_id_i <- as.numeric(readline("Bitte zu zutreffende ID ZWISCHEN 1 UND 1 WÄHLEN! ->   "))
    }
    
    cat("\nWelche Masse wurde eingewogen?\n")
    extraction_mass_i <- as.numeric(readline("Bitte die eingewogene Masse in g angeben ->   "))
    while (is.na(extraction_mass_i)){
      cat("\nKeine gültige Zahl gefunden! Gegebenenfalls , durch . ersetzen.\n")
      extraction_mass_i <- as.numeric(readline("Bitte die eingewogene Masse in g angeben! ->   "))
    }
    leaching_id_i <- 0}
  
  if (sample_type_ii==2) {
    sample_type_i <- "Leaching"
    cat("\nWelche Leachingmethode wurde verwendet?
    [1] Solvent= MQ - Zeit= 9 d - T= RT \n")
    leaching_id_i <- as.numeric(readline("Bitte zutreffende ID in die Console eingeben und bestätigen ->   "))
    while (!any(leaching_id_i==(c(1:1))) | is.na(leaching_id_i)){
      cat("\nDie angegebene ID steht nicht zur Auswahl!\n")
      leaching_id_i <- as.numeric(readline("Bitte zu zutreffende ID ZWISCHEN 1 UND 1 WÄHLEN! ->   "))
    }
    cat("\nWurde die Masse oder die Oberfläche der Probe erhoben?
        [1] Masse
        [2] Oberfläche\n")
    leaching_mass_area_i <- as.numeric(readline("Bitte die zutreffende ID in die Console eingeben ->   "))
    while (!any(leaching_mass_area_i==(c(1:2))) | is.na(leaching_mass_area_i)){
      cat("\nDie angegebene ID steht nicht zur Auswahl!\n")
      leaching_mass_area_i <- as.numeric(readline("Bitte zu zutreffende ID ZWISCHEN 1 UND 2 WÄHLEN! ->   "))
    }
    if (leaching_mass_area_i==1){
      leaching_mass_i <- as.numeric(readline("Bitte die eingewogene Masse in g angeben ->   "))
      while (is.na(leaching_mass_i)){
        cat("\nKeine gültige Zahl gefunden! Gegebenenfalls , durch . ersetzen.\n")
        leaching_mass_i <- as.numeric(readline("Bitte die eingewogene Masse in g angeben! ->   "))
      }
      leaching_area_i <- 0
      }
    if (leaching_mass_area_i==2){
      leaching_area_i <- as.numeric(readline("Bitte die Oberfläche der Probe angeben ->   "))
      while (is.na(leaching_area_i)){
        cat("\nKeine gültige Zahl gefunden! Gegebenenfalls , durch . ersetzen.\n")
        leaching_area_i <- as.numeric(readline("Bitte die Oberfläche der Probe angeben! ->   "))
      }
      leaching_mass_i <- 0
      }
    extraction_id_i <- 0
    }
  
  cat("\nWoher stammt die Probe?\n")
  location_i <- readline("Bitte Ortsbezeichnung in die Console eingeben und bestätigen ->   ")
  while (location_i=="" | !is.na(as.integer(location_i))) {
    cat("\nBitte Ort angeben! Stammt die Probe vielleicht aus dem Labor?\n")
    location_i <- readline("Bitte einen anderen Probenname angeben ->   ")
  }
  
  cat("\nWurde die Probe angereichert?
    Wenn ja, bitte den Anreicherungsfaktor in die Console eingeben
    Wenn nein, bitte 0 in die Console eingeben\n")
  enrichment_factor_i <- as.numeric(readline("->   "))
  while(is.na(enrichment_factor_i)){ 
    cat("\nBitte den Anreicherungsfaktor angeben!
        Wenn die Probe nicht angereichert wurde, dann 0 eingeben!\n")
    enrichment_factor_i <- as.numeric(readline("->   "))
    }
  if (enrichment_factor_i>0) {enrichment <- TRUE
  } else {enrichment <- FALSE}
  
  cat("\nWann wurde die Probe erstellt?\n")
  date_i <- readline("Bitte Datum im Format YYYY-MM-DD angeben und bestätigen ->   ")
  while (is.na(as.Date(stringr::str_extract(date_i, "\\d{4}-\\d{2}-\\d{2}"), format="%Y-%m-%d")) |
         as.Date(stringr::str_extract(date_i, "\\d{4}-\\d{2}-\\d{2}"), format="%Y-%m-%d") > Sys.Date()){
    if (!is.na(as.Date(stringr::str_extract(date_i, "\\d{4}-\\d{2}-\\d{2}"), format="%Y-%m-%d")) &&
      as.Date(stringr::str_extract(date_i, "\\d{4}-\\d{2}-\\d{2}"), format="%Y-%m-%d") > Sys.Date()) {cat("Bist du Marty McFly?")}
    cat("\nKeine gültige Eingabe!\n")
    date_i <- readline("Bitte Datum im Format YYYY-MM-DD angeben und bestätigen ->   ")
  }
    
  cat("\nWer ist der Ansprechpartner für eventuelle Rückfragen zu dieser Probe?\n")
  contact_i <- readline("Bitte den Namen eingeben und bestätigen ->   ")
  while (contact_i=="" | !is.na(as.integer(contact_i))) {
    cat("\nWer ist der Ansprechpartner für eventuelle Rückfragen zu dieser Probe?\n")
    contact_i <- readline("Bitte einen Ansprechpartner angeben! ->   ")
  }
  
  cat("\nKann die Probe einem Projekt zugeordnet werden?
    Wenn nein, dann ohne Texteingabe bestätigen\n")
  project_i <- readline("Bitte den Projektnamen eingeben und bestätigen ->   ")
  
  cat("\nSoll die Json direkt in die Label_DB importiert werden?
      [1]  Ja
      [2]  Nein")
  import_i <- as.numeric(readline("Bitte zutreffende ID in die Console eingeben und bestätigen ->   "))
  while (!any(import_i==(c(1:2))) | is.na(import_i)){
    cat("\nDie angegebene ID steht nicht zur Auswahl!\n")
    import_i <- as.numeric(readline("Bitte zu zutreffende ID ZWISCHEN 1 UND 2 WÄHLEN! ->   "))
  }
  
  
  cat("\nProcessing, please wait!\n")
  
  
  # load Spectra_DB and Label_DB
  DB_Spectra <- DBI::dbConnect(RSQLite::SQLite(), paste0(db_path,"/MS2_db_v7.db"))
  DB_Lable <- DBI::dbConnect(RSQLite::SQLite(), paste0(db_path,"/MS2_db_v5.db"))
  #DB_Spectra <- DBI::dbConnect(RSQLite::SQLite(), "F:\\Homeoffice\\202102089\\MS2_db_v7.db")
  #DB_Lable <- DBI::dbConnect(RSQLite::SQLite(), "F:\\Homeoffice\\202102089\\MS2_db_v5.db")
  
  
  # looking for other Json files and pick the max IDs from the latest json file
  old_json_files <- list.files(json_path,pattern=".json$",full.names = TRUE,recursive = TRUE)
  #old_json_files <- list.files(path="F:\\Homeoffice\\202102089",pattern=".json$",full.names = TRUE,recursive = TRUE)
  if (!length(old_json_files)==0){
    dirs <- dirname(old_json_files)
    lastfiles <- tapply(old_json_files,dirs,function(v) v[which.max(file.mtime(v))])
    latest_json_file <- jsonlite::fromJSON(lastfiles)
    
    experiment_id <- max(latest_json_file$experiment$experiment_id)
    fragment_id <- max(latest_json_file$fragment$fragment_id)
    sample_id<- max(latest_json_file$sample$sample_id)+1
    
  }  else {
    experiment_id <- 0
    fragment_id <- 0
    sample_id<- 1
  }     
  
  ## Prozessierung der Daten aus der App ##  
  
  #browser()
  
  for (i in 1:nrow(grouped_i)){
    
    # show Progress in 10% steps
    if (stringr::str_detect(as.character(as.integer(i/nrow(grouped_i)*100)),"0")){
      print(paste(as.integer(i/nrow(grouped_i)*100),"%"))}
    
    # get mz, rt and peaklist/ Data for experiment column
    # find sample/ peaklist with highest intensity
    col_high_int <- which.max(grouped_i[i,grep("Int_",colnames(grouped_i))[samples_i]])  
    PeakID_high_int <- grouped_i[i,sprintf("PeakID_%i",col_high_int)]  
    sample_high_int <- col_high_int
    
    Zeile <- peaklist_i[[sample_high_int]]$peak_id_all==PeakID_high_int
    mz <- round(peaklist_i[[sample_high_int]][Zeile,]$mz,4)
    RT <- round(peaklist_i[[sample_high_int]][Zeile,]$RT/60,2)
    
    # Get MS2/ Data for fragment column
    if (peaklist_i[[sample_high_int]]$MS2scan[Zeile]==0)  next
    ms2 <-xcms::getMsnScan(datenlist_i[[sample_high_int]], peaklist_i[[sample_high_int]]$MS2scan[Zeile])
    ms2 <-as.data.frame(ms2, row.names = c("mz","intensity"))
    ms2 <-ms2[ms2$intensity> max(ms2$intensity)*0.02,]
    ms2 <-ms2[ms2$mz < mz+1,]
    
    # get annotated compounds/ Data for compound column
    if (exists("annotationTable_i")){
      alignmentID <- grouped_i[i,"alignmentID"]
      if (any(annotationTable_i$alignmentID==alignmentID)){
        #browser()
        name_i <-annotationTable_i[annotationTable_i$alignmentID==alignmentID,"name"]
        CAS_i <- annotationTable_i[annotationTable_i$alignmentID==alignmentID,"CAS"]
        formula_i <- annotationTable_i[annotationTable_i$alignmentID==alignmentID,"formula"]
        SMILES_i <- annotationTable_i[annotationTable_i$alignmentID==alignmentID,"SMILES"]
        adduct_i <- annotationTable_i[annotationTable_i$alignmentID==alignmentID,"adduct"]
        isotope_i <- annotationTable_i[annotationTable_i$alignmentID==alignmentID,"isotope"]
        Spectra_DB_score_i <- annotationTable_i[annotationTable_i$alignmentID==alignmentID,"score"]
      } else {
        name_i <- NA
        CAS_i <- NA
        formula_i <- NA
        SMILES_i <- NA
        adduct_i <- NA
        isotope_i <- NA
        Spectra_DB_score_i <- NA
      }}
    
    
    ## Filling Data for every Feature in lists ##
    
    # Compound
    if (is.na(name_i)) {compound_id <- 0
    
    } else {
      compound_match_Spectra <- DBI::dbReadTable(DB_Spectra, "compound") %>% dplyr::filter(name==name_i[[1]])
      compound_match_Lable <- DBI::dbReadTable(DB_Lable, "compound") %>% dplyr::filter(name==name_i[[1]])
      
      if (nrow(compound_match_Spectra)!=0 && nrow(compound_match_Lable)==0){
        compound_id <- compound_match_Spectra[1]
        compound <- rbind(compound,compound_match_Spectra)
        
      } else {compound_id=compound_match_Spectra[1]}
    }
    
    # Experiment and expGroupExp
    experiment_id = experiment_id+1
    #sample_id_experiment <- c(sample_id:(sample_id+repikat_i-1))[col_high_int]
    sample_id_experiment <- c(sample_id:(sample_id+length(samples_i)-1))[col_high_int]
    new_experiment <- data.frame(experiment_id = experiment_id, 
                                 compound_id=compound_id, 
                                 parameter_id=parameter_id_i, 
                                 sample_id=sample_id_experiment, 
                                 mz=mz,
                                 rt=RT,
                                 time_added=Sys.Date(), 
                                 adduct=adduct_i, 
                                 isotope=isotope_i,
                                 Spectra_DB_score=Spectra_DB_score_i[[1]])
    experiment <- rbind(experiment,new_experiment)
    
    new_expGroupExp <- data.frame(experimentGroup_id=experimentGroup_id_i,
                                  experiment_id=experiment_id                                  )
    expGroupExp <- rbind(expGroupExp,new_expGroupExp)
    
    # Fragment
    for (j in 1:nrow(ms2)){
      fragment_id=fragment_id+1
      new_fragment <- data.frame(fragment_id=fragment_id,
                                 mz=ms2$mz[j],
                                 int=ms2$intensity[j],
                                 experiment_id = experiment_id)
      fragment <- rbind(fragment,new_fragment)
    }
    
  }
  
  
  ## Filling Data for the sample in the list ##
  
  # sample <- data.frame(sample_id=c(sample_id:(sample_id+repikat_i-1)), 
  sample <- data.frame(sample_id=c(sample_id:(sample_id+length(samples_i)-1)),
                       samplename=samplename_i, 
                       sample_type=sample_type_i, 
                       location=location_i, 
                       enrichment=enrichment, 
                       enrichment_factor=enrichment_factor_i, 
                       data_name=stringr::str_match(basename(samplelist_i$File[samples_i]), "(.*)\\.mzXML")[,2], 
                       data_location=samplelist_i$File[samples_i], 
                       date=date_i, 
                       contact=contact_i, 
                       project=project_i, 
                       leaching_id=leaching_id_i, 
                       extraction_id=extraction_id_i,
                       stringsAsFactors=FALSE)
  
  # sampLeachSamp and samplExtractSamp
  if (leaching_id_i!=0){
    sampLeachSamp <- data.frame(sample_id=c(sample_id:(sample_id+length(samples_i)-1)), 
                                leaching_id=leaching_id_i)
    leaching$sample_mass <- leaching_mass_i
    leaching$sample_area <- leaching_area_i
    
  } else {samplExtractSamp <- data.frame(sample_id=c(sample_id:(sample_id+length(samples_i)-1)), 
                                         extraction_id=extraction_id_i)
  extraction$sample_mass <- extraction_mass_i
  }
  
  
  ## Save data as Json und disconnect DBs ##
  
  DB_data <- list(compGroupComp,compound,compoundGroup,experiment,experimentGroup,expGroupExp,extraction,fragment,
                  leaching,parameter,sample,sampLeachSamp,samplExtractSamp)
  names(DB_data) <- c("compGroupComp","compound","compoundGroup","experiment","experimentGroup","expGroupExp","extraction","fragment",
                      "leaching","parameter","sample","sampLeachSamp","samplExtractSamp")
  
  cat("\nJSON Finshed!\n
      Json-File ",samplename_i," saved\n")
  
  DB_data_json <- jsonlite::toJSON(DB_data)
  write(DB_data_json,file=paste0(json_path,"/",samplename_i,".json"))
  #write(DB_data_json,file=paste0("F:\\Homeoffice\\202102089\\",samplename_i,".json"))
  
  DBI::dbDisconnect(DB_Lable)
  DBI::dbDisconnect(DB_Spectra)
  
  
  # If desired, the created Json is imported directly into the DB 
  if (import_i==1) {insert_json_data()}
  
}



#' insert_json_data
#'
#' @param json_path path to json-files on Z
#' @param db_path path to both DBs on Server (Spectra_DB = db_v7 and Label_DB = db_v5)
#'
#' @return insert json-files to Label_DB (db_v5)
#' @export
#' @import DBI
#' @import RSQLite
#'
#' @examples
insert_json_data <- function(json_path="~/Data/Toni/Lable_Datenbank/Json_Files",
                             db_path="~/sqlite_local"){
  
  stopifnot(file.exists("~/sqlite_local/MS2_db_v5.db"))
  #stopifnot(file.exists("~/Data/Toni/Lable_Datenbank/Test_DB/MS2_db_v5.db"))
  
  # load Label_DB and Json-files
  
  DB_Lable <- DBI::dbConnect(RSQLite::SQLite(), paste0(db_path,"/MS2_db_v5.db"))
  #DB_Lable <- DBI::dbConnect(RSQLite::SQLite(), "~/Data/Toni/Lable_Datenbank/Test_DB/MS2_db_v5.db")
  
  all_json_files <- list.files(path=json_path,pattern=".json$",full.names = TRUE,recursive = TRUE)
  
  #stopifnot(file.exists("F:\\Homeoffice\\202102089\\MS2_db_v5.db"))
  #DB_Lable <- DBI::dbConnect(RSQLite::SQLite(), "F:\\Homeoffice\\202102089\\MS2_db_v5.db")
  #all_json_files <- list.files(path="F:\\Homeoffice\\202102089",pattern=".json$",full.names = TRUE,recursive = TRUE)
  
  #insertet_json.txt containes json-files, already been inserted to Label_DB
  
  if(file.exists(paste0(json_path,"/insertet_json.txt"))){
    inserted_json <- read.table(paste0(json_path,"/insertet_json.txt"),stringsAsFactors = F)[,1]
  } else {
    inserted_json <- vector()
  }
  
  #  if(file.exists("F:\\Homeoffice\\202102089\\insertet_json.txt")){
  #    inserted_json <- read.table("F:\\Homeoffice\\202102089\\insertet_json.txt",stringsAsFactors = F)[,1]
  #  } else {
  #    inserted_json <- vector()
  #  }
  
  # look for new json-files
  if(length(inserted_json)!=0) {
    new_json_files <- vector()
    new_json_files <- setdiff(all_json_files,inserted_json)
    dirs <- dirname(new_json_files)
  } else {
    new_json_files <- all_json_files
    dirs <- dirname(new_json_files)
  }
  
  if(length(new_json_files)==0){ 
    stop("Es liegen keine neuen JSON-Files vor.\n
         Die Datenbank ist auf dem neuesten Stand")}
  
  good_json <- vector()
  bad_json <- vector()
  
  
  # loop to insert every single json
  for (i in 1:length(new_json_files)){
    
    ordered_files <- tapply(new_json_files,dirs,function(v) v[order(file.mtime(v))])[[1]]
    oldest_file_json_file <- jsonlite::fromJSON(ordered_files[i])
    
    DBI::dbBegin(DB_Lable)
    
    tryCatch({
      
      #new feature-specific data
      DBI::dbWriteTable(DB_Lable, "fragment", oldest_file_json_file$fragment, append= T, overwrite=F)
      DBI::dbWriteTable(DB_Lable, "experiment", oldest_file_json_file$experiment, append= T, overwrite=F)
      DBI::dbWriteTable(DB_Lable, "sample", oldest_file_json_file$sample, append= T, overwrite=F)
      
      
      # new relationship-data
      if (length(oldest_file_json_file$sampLeachSamp)!=0){
        DBI::dbWriteTable(DB_Lable, "sampLeachSamp", oldest_file_json_file$sampLeachSamp, append= T, overwrite=F)}
      
      if (length(oldest_file_json_file$samplExtractSamp)!=0){
        DBI::dbWriteTable(DB_Lable, "samplExtractSamp", oldest_file_json_file$samplExtractSamp, append= T, overwrite=F)}
      
      if (length(oldest_file_json_file$compGroupComp)!=0){
        DBI::dbWriteTable(DB_Lable, "compGroupComp", oldest_file_json_file$compGroupComp, append= T, overwrite=F)}
      
      if (length(oldest_file_json_file$expGroupExp)!=0){
        DBI::dbWriteTable(DB_Lable, "expGroupExp", oldest_file_json_file$expGroupExp, append= T, overwrite=F)}
      
      
      # new catalogue-data
      # wenn keine annotationtable vorhanden war ist compound NULL wodurch keine neuen Werte eingetragen werden
      if  (!is.null(nrow(oldest_file_json_file$compound))){
        compound_DB <- DBI::dbReadTable(DB_Lable, "compound")
        new_compound <-  oldest_file_json_file$compound[!(oldest_file_json_file$compound$compound_id %in% compound_DB$compound_id),]
        if (nrow(new_compound)>0) {
          DBI::dbWriteTable(DB_Lable, "compound", new_compound, append= T, overwrite=F)}
      }
      
      compoundGroup_DB <- DBI::dbReadTable(DB_Lable, "compoundGroup")
      new_compoundGroup <-  oldest_file_json_file$compoundGroup[!(oldest_file_json_file$compoundGroup$compoundGroup_id %in% compoundGroup_DB$compoundGroup_id),]
      if (nrow(new_compoundGroup)>0) {
        DBI::dbWriteTable(DB_Lable, "compoundGroup", new_compoundGroup, append= T, overwrite=F)}
      
      extraction_DB <- DBI::dbReadTable(DB_Lable, "extraction")
      new_extraction <-  oldest_file_json_file$extraction[!(oldest_file_json_file$extraction$extraction_id %in% extraction_DB$extraction_id),]
      if (nrow(new_extraction)>0) {
        DBI::dbWriteTable(DB_Lable, "extraction", new_extraction, append= T, overwrite=F)}
      
      leaching_DB <- DBI::dbReadTable(DB_Lable, "leaching")
      new_leaching <-  oldest_file_json_file$leaching[!(oldest_file_json_file$leaching$leaching_id %in% leaching_DB$leaching_id),]
      if (nrow(new_leaching)>0) {
        DBI::dbWriteTable(DB_Lable, "leaching", new_leaching, append= T, overwrite=F)}
      
      parameter_DB <- DBI::dbReadTable(DB_Lable, "parameter")
      new_parameter <-  oldest_file_json_file$parameter[!(oldest_file_json_file$parameter$parameter_id %in% parameter_DB$parameter_id),]
      if (nrow(new_parameter)>0) {
        DBI::dbWriteTable(DB_Lable, "parameter", new_parameter, append= T, overwrite=F)}
      
      experimentGroup_DB <- DBI::dbReadTable(DB_Lable, "experimentGroup")
      new_experimentGroup <-  oldest_file_json_file$experimentGroup[!(oldest_file_json_file$experimentGroup$experimentGroup_id %in% experimentGroup_DB$experimentGroup_id),]
      if (nrow(new_experimentGroup)>0) {
        DBI::dbWriteTable(DB_Lable, "experimentGroup", new_experimentGroup, append= T, overwrite=F)}
      
      # write processed json-file in insert txt and comit insertion of actual json-file
      #inserted_json[(nrow(inserted_json)+1),] <- ordered_files[i]
      good_json <- c(good_json,ordered_files[i])
      DBI::dbCommit(DB_Lable)
      
      # error and end of trycatch
      
    },
    # If an error occurs, the rollback command resets the database to the state before
    # this file was inserted.
    error=function(x) {
      DBI::dbRollback(DB_Lable)
      bad_json <- c(bad_json,ordered_files[i])
      print(paste("Fehler mit dem Json_File", ordered_files[i]))
      
    }) 
  }
  # update inserted_json.txt and disconnect database
  #write.table(inserted_json,"F:\\Homeoffice\\202102089\\insertet_json.txt",row.names = F,col.names = F)
  if (length(good_json)!=0){
    inserted_json <- c(inserted_json,good_json)
    write.table(inserted_json, file= paste0(json_path,"/insertet_json.txt"),row.names = F,col.names = F)
    print(paste("Die Json-Files",good_json, "wurden in die Label_DB importiert."))
  } else {
    print("Es wurden keine neuen Daten in die Label_DB importiert.")
  }
  DBI::dbDisconnect(DB_Lable)
  
}                                  




