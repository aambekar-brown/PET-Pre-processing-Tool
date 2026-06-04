
library(shiny)
library(fslr)
library(neurobase)
library(oro.nifti)
library(grid)
library(extrantsr)
library(neuroCombat)
library(RAVEL)
library(ANTsR)

library(shinyFiles)
library(reticulate)
library(rgl)
library(png)
library(readxl)

source("ReportGenerator.R")

#install.packages(c("shinyFiles", "reticulate", "rgl", "png", "readxl"))


use_python("/usr/bin/python3", required = TRUE)

# Define the user interface
ui <- fluidPage(
  tags$head(
    tags$style(HTML("
    /* ===== Increase overall UI font size ===== */
    body { font-size: 20px; }

    /* Inputs + labels */
    .shiny-input-container { font-size: 20px; }
    .shiny-input-container label { font-size: 20px; }
    .form-control { font-size: 20px; }

    /* Buttons */
    .btn { font-size: 20px; }

    /* Headings (optional) */
    h4 { font-size: 25px; }

    .preprocessing-steps {
      min-height:260px
      height: auto;
    }
    .preprocessing-steps .checkbox {
      margin-bottom: 10px;
    }
    .console-container {
      min-height:260px
      height: auto;
    }
  "))
  ),
  titlePanel(
    title = tags$span(tags$b("TRAECR"), ": PET/MRI Pre-processing Tool"),
    windowTitle = "PET Pre-processing Tool"
  ),
  
  wellPanel(
    div(
      fluidRow(
        column(12,
               actionButton("artifact_viewer", "NIfTI Image Properties QC Dashboard (MRI/PET)", class = "btn-primary", style = "width:100%; margin-bottom: 20px; font-weight: bold;")
        )
      ),
      fluidRow(
        column(3, strong("PET Scan Input Path")),
        column(6, textInput("pet_input", label = NULL, placeholder = "Path to PET scans")),
        column(3, shinyFilesButton("btn_pet", "Select", title = "Select folder", multiple = TRUE))
      ),
      fluidRow(
        column(3, strong("MRI Scans Input Path")),
        column(6, textInput("mri_input", label = NULL, placeholder = "Path to MRI scans")),
        column(3, shinyFilesButton("btn_mri", "Select", title = "Select folder", multiple = TRUE))
      ),
      fluidRow(
        column(3, strong("CSV File Input Path")),
        column(6, textInput("csv_input", label = NULL, placeholder = "Path to a CSV or Excel file")),
        column(3, shinyFilesButton("btn_csv", "Select", title = "Select CSV file", multiple = FALSE))
      ),
      fluidRow(
        column(3, strong("Template Selection")),
        column(6, selectInput("template_selection", label = NULL, choices = c("MNI152_T1_1mm.nii", "MNI152_T1_2mm.nii"))),
      ),
      fluidRow(
        column(3, strong("Control Region Mask Input Path")),
        column(6, textInput("control_mask_input", label = NULL, placeholder = "Path to control region mask")),
        column(3, shinyFilesButton("btn_control_mask", "Select", title = "Select control region mask", multiple = FALSE))
      ),
      fluidRow(
        column(3, strong("Brain Mask Input Path")),
        column(6, textInput("brain_mask_input", label = NULL, placeholder = "Path to brain region mask")),
        column(3, shinyFilesButton("btn_brain_mask", "Select", title = "Select a brain mask", multiple = FALSE))
      )
    )
  ),
  
  fluidRow(
    column(4,
           wellPanel(
             div(class = "preprocessing-steps",
                 radioButtons("preprocessing_steps",
                              h4(strong("Preprocessing Type")),
                              choices = list("Brain Extraction" = "brain_extraction",
                                             "Template Registration" = "mni_registration",
                                             "MRI PET Co-registration" = "mri_coregistration",
                                             "COMBAT Harmonization" = "combat_harmonization",
                                             "RAVEL Normalization" = "ravel_normalization")
                 ),
                 
                 conditionalPanel(
                   condition = "input.preprocessing_steps == 'mni_registration' || input.preprocessing_steps == 'mri_coregistration'",
                   checkboxInput("do_parcellation", "(Optional) Parcellation", value = FALSE)
                 ),
                 actionButton("preprocess", "Preprocess", class = "btn-primary", style = "font-weight: bold;")
             )
           )
    ),
    column(8,
           wellPanel(
             div(class = "console-container",
                 h4(strong("Output Console")),
                 div(style = "height: 200px; overflow-y: auto; background-color: #f7f7f9; border-radius: 5px; padding: 5px;",
                     textOutput("console", container = tags$pre))
             )
           )
    )
  ),
  fluidRow(
    column(12,
           wellPanel(
             div(class = "visualization-container",
                 h4(strong("Output Visualization")),
                 selectInput("selected_output", "Select Output:", choices = NULL),
                 div(style = "display: flex; align-items: center; gap: 10px;",
                     actionButton("prev_slice", "Previous Slice", style = "flex: 1;"),
                     sliderInput("slice_number", "Slice Number:", min = 1, max = 1, value = 1, width = "50%"),
                     actionButton("next_slice", "Next Slice", style = "flex: 1;"),
                     actionButton("play_slices", "Play", style = "flex: 1;"),
                     actionButton("stop_slices", "Pause", style = "flex: 1;")
                 ),
                 plotOutput("brain_plot", width = "100%", height = "1000px")
             )
           )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Create a reactive value to store the console messages
  consoleMessages <- reactiveVal("")
  
  
  # -----------------------
  # Optional parcellation
  # -----------------------
  .find_default_atlas <- function() {
    tpl_dir <- file.path(getwd(), "templates")
    if (!dir.exists(tpl_dir)) return(NA_character_)
    cand <- list.files(tpl_dir, pattern = "(AAL|atlas|parcell|label).*\\.nii(\\.gz)?$", ignore.case = TRUE, full.names = TRUE)
    if (length(cand) == 0) return(NA_character_)
    cand[[1]]
  }
  
  run_parcellation <- function(
    target_img_path,
    output_dir,
    prefix = "parcellation",
    # NEW (optional) inputs for proper MNI->PET workflow
    mri_to_mni_mat = NULL,     # MRI->MNI (from your mri.mat when ref=template)
    pet_to_mri_mat = NULL,     # PET->MRI (from flirt omat when ref=MRI)
    pet_ref_path   = NULL,     # ORIGINAL PET path (defines "PET space")
    mri_ref_path   = NULL      # optional MRI ref if PET not available
  ) {
    atlas_path <- .find_default_atlas()
    if (!is.character(atlas_path) || !nzchar(atlas_path) || is.na(atlas_path) || !file.exists(atlas_path)) {
      msg <- paste0(
        "Parcellation skipped: no atlas NIfTI found under ",
        file.path(getwd(), "templates"),
        " (expected AAL*.nii(.gz) / *atlas*.nii(.gz) / *label*.nii(.gz)).\n"
      )
      consoleMessages(paste0(consoleMessages(), msg))
      return(list(performed = FALSE, note = msg))
    }
    
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    # -----------------------------
    # Decide what "PET space" is
    # -----------------------------
    ref_path <- NA_character_
    if (is.character(pet_ref_path) && nzchar(pet_ref_path) && file.exists(pet_ref_path)) {
      ref_path <- pet_ref_path
    } else if (is.character(mri_ref_path) && nzchar(mri_ref_path) && file.exists(mri_ref_path)) {
      ref_path <- mri_ref_path
    } else if (is.character(target_img_path) && nzchar(target_img_path) && file.exists(target_img_path)) {
      # fallback: whatever target was passed
      ref_path <- target_img_path
    }
    
    if (!is.character(ref_path) || is.na(ref_path) || !file.exists(ref_path)) {
      msg <- "Parcellation skipped: no valid reference image path found for resampling.\n"
      consoleMessages(paste0(consoleMessages(), msg))
      return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
    }
    
    # -----------------------------
    # If we have transforms:
    #   MRI->MNI and PET->MRI
    # then compute:
    #   MNI->MRI (inverse)
    #   MRI->PET (inverse of PET->MRI)
    #   MNI->PET = (MRI->PET) o (MNI->MRI)
    # and resample atlas into PET space with NN.
    # -----------------------------
    atlas_in_ref <- NULL
    atlas_in_ref_path <- file.path(output_dir, paste0(prefix, "_atlas_in_ref.nii.gz"))
    
    have_mri_to_mni <- isTRUE(!is.null(mri_to_mni_mat) && length(mri_to_mni_mat) == 1 &&
                                is.character(mri_to_mni_mat) && nzchar(mri_to_mni_mat) &&
                                file.exists(mri_to_mni_mat))
    
    have_pet_to_mri <- isTRUE(!is.null(pet_to_mri_mat) && length(pet_to_mri_mat) == 1 &&
                                is.character(pet_to_mri_mat) && nzchar(pet_to_mri_mat) &&
                                file.exists(pet_to_mri_mat))
    
    if (have_mri_to_mni && have_pet_to_mri) {
      # Build needed mats
      mni_to_mri_mat <- file.path(output_dir, paste0(prefix, "_mni2mri.mat"))
      mri_to_pet_mat <- file.path(output_dir, paste0(prefix, "_mri2pet.mat"))
      mni_to_pet_mat <- file.path(output_dir, paste0(prefix, "_mni2pet.mat"))
      
      # Invert MRI->MNI to get MNI->MRI
      system2("convert_xfm", args = c("-inverse", "-omat", mni_to_mri_mat, mri_to_mni_mat), stdout = TRUE, stderr = TRUE)
      
      # Invert PET->MRI to get MRI->PET
      system2("convert_xfm", args = c("-inverse", "-omat", mri_to_pet_mat, pet_to_mri_mat), stdout = TRUE, stderr = TRUE)
      
      # Compose: MNI->PET = (MRI->PET) o (MNI->MRI)
      # For FSL: -concat A B gives A*B (apply B then A)
      system2("convert_xfm", args = c("-omat", mni_to_pet_mat, "-concat", mri_to_pet_mat, mni_to_mri_mat),
              stdout = TRUE, stderr = TRUE)
      
      if (!file.exists(mni_to_pet_mat)) {
        msg <- "Parcellation skipped: failed to create composed MNI->PET transform matrix.\n"
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
      # Resample atlas into PET space (or chosen ref space) using NN
      # NOTE: ref_path is PET if provided, otherwise MRI/target fallback
      system2("flirt", args = c(
        "-in", atlas_path,
        "-ref", ref_path,
        "-applyxfm",
        "-init", mni_to_pet_mat,
        "-interp", "nearestneighbour",
        "-out", atlas_in_ref_path
      ), stdout = TRUE, stderr = TRUE)
      
      if (!file.exists(atlas_in_ref_path)) {
        msg <- "Parcellation skipped: failed to resample atlas into PET/reference space.\n"
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
      atlas_in_ref <- tryCatch(readNIfTI(atlas_in_ref_path, reorient = FALSE), error = function(e) NULL)
      if (is.null(atlas_in_ref)) {
        msg <- "Parcellation skipped: could not read resampled atlas in PET/reference space.\n"
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
    } else if (have_mri_to_mni && !have_pet_to_mri) {
      
      # MRI-only: resample atlas from MNI into MRI/reference space using inverse(mri_to_mni_mat)
      mni_to_mri_mat <- file.path(output_dir, paste0(prefix, "_mni2mri.mat"))
      
      system2("convert_xfm",
              args = c("-inverse", "-omat", mni_to_mri_mat, mri_to_mni_mat),
              stdout = TRUE, stderr = TRUE)
      
      if (!file.exists(mni_to_mri_mat)) {
        msg <- "Parcellation skipped: failed to invert MRI->MNI matrix (cannot create MNI->MRI).\n"
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
      system2("flirt", args = c(
        "-in", atlas_path,
        "-ref", ref_path,
        "-applyxfm",
        "-init", mni_to_mri_mat,
        "-interp", "nearestneighbour",
        "-out", atlas_in_ref_path
      ), stdout = TRUE, stderr = TRUE)
      
      if (!file.exists(atlas_in_ref_path)) {
        msg <- "Parcellation skipped: failed to resample atlas into MRI/reference space.\n"
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
      atlas_in_ref <- tryCatch(readNIfTI(atlas_in_ref_path, reorient = FALSE), error = function(e) NULL)
      if (is.null(atlas_in_ref)) {
        msg <- "Parcellation skipped: could not read resampled atlas in MRI/reference space.\n"
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
    }
    else {
      # No transform info -> we cannot do proper MNI->PET workflow.
      # Keep behavior safe: resample atlas to ref grid ONLY IF it already matches
      # (otherwise skip). This prevents wrong ROI outputs.
      ref_img <- tryCatch(readNIfTI(ref_path, reorient = FALSE), error = function(e) NULL)
      atlas_img <- tryCatch(readNIfTI(atlas_path, reorient = FALSE), error = function(e) NULL)
      
      if (is.null(ref_img) || is.null(atlas_img)) {
        msg <- "Parcellation skipped: could not read reference image or atlas image.\n"
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
      if (!all(dim(ref_img)[1:3] == dim(atlas_img)[1:3])) {
        msg <- paste0(
          "Parcellation skipped: transforms not provided and atlas grid != reference grid. Ref dims=",
          paste(dim(ref_img)[1:3], collapse="x"),
          ", atlas dims=", paste(dim(atlas_img)[1:3], collapse="x"),
          ".\n"
        )
        consoleMessages(paste0(consoleMessages(), msg))
        return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
      }
      
      # If same grid, just use the atlas directly
      atlas_in_ref <- atlas_img
      # still write it out consistently
      tryCatch(writeNIfTI(atlas_in_ref, filename = sub("\\.nii\\.gz$", "", atlas_in_ref_path)),
               error = function(e) NULL)
    }
    
    # -----------------------------
    # Compute ROI means on PET (ref_path)
    # -----------------------------
    ref_img <- tryCatch(readNIfTI(ref_path, reorient = FALSE), error = function(e) NULL)
    if (is.null(ref_img)) {
      msg <- "Parcellation skipped: could not read PET/reference image for ROI stats.\n"
      consoleMessages(paste0(consoleMessages(), msg))
      return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
    }
    
    ref_vals  <- as.numeric(ref_img)
    atlas_vals <- as.integer(round(as.numeric(atlas_in_ref)))
    atlas_vals[!is.finite(atlas_vals)] <- 0
    
    ok <- is.finite(ref_vals) & is.finite(atlas_vals) & atlas_vals > 0
    if (!any(ok)) {
      msg <- "Parcellation skipped: no labeled voxels (atlas > 0) overlapped finite voxels in PET/reference.\n"
      consoleMessages(paste0(consoleMessages(), msg))
      return(list(performed = FALSE, note = msg, atlas_path = atlas_path))
    }
    
    lbls <- sort(unique(atlas_vals[ok]))
    means <- vapply(lbls, function(L) mean(ref_vals[ok & atlas_vals == L], na.rm = TRUE), numeric(1))
    nvox  <- vapply(lbls, function(L) sum(ok & atlas_vals == L), integer(1))
    
    parc_df <- data.frame(Label = lbls, Nvox = nvox, Mean = means, stringsAsFactors = FALSE)
    
    out_csv <- file.path(output_dir, paste0(prefix, "_roi_means.csv"))
    tryCatch(write.csv(parc_df, out_csv, row.names = FALSE), error = function(e) NULL)
    
    # Save atlas labels map (already in PET/reference grid)
    atlas_labels_base <- file.path(output_dir, paste0(prefix, "_atlas_labels_in_ref"))
    tryCatch({
      atlas_out <- atlas_in_ref
      atlas_out[] <- atlas_vals
      writeNIfTI(atlas_out, filename = atlas_labels_base)  # writes .nii.gz
    }, error = function(e) {
      msg2 <- paste0("Warning: could not write atlas label-map: ", e$message, "\n")
      consoleMessages(paste0(consoleMessages(), msg2))
    })
    
    # Save one binary ROI mask per label (in PET/reference grid)
    roi_dir <- file.path(output_dir, paste0(prefix, "_rois_in_ref"))
    dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)
    
    for (L in lbls) {
      roi_base <- file.path(roi_dir, paste0(prefix, "_roi_", L))
      tryCatch({
        roi_img <- atlas_in_ref
        roi_img[] <- as.integer(atlas_vals == L)
        writeNIfTI(roi_img, filename = roi_base)  # writes .nii.gz
      }, error = function(e) {
        msg2 <- paste0("Warning: could not write ROI mask for label ", L, ": ", e$message, "\n")
        consoleMessages(paste0(consoleMessages(), msg2))
      })
    }
    
    msg <- paste0(
      "Parcellation completed in PET/reference space using atlas: ", basename(atlas_path),
      "\n  Resampled atlas: ", atlas_in_ref_path,
      "\n  ROI means: ", out_csv,
      "\n  ROI masks dir: ", roi_dir, "\n"
    )
    consoleMessages(paste0(consoleMessages(), msg))
    
    return(list(
      performed = TRUE,
      atlas_path = atlas_path,
      ref_path = ref_path,
      resampled_atlas_path = atlas_in_ref_path,
      csv_path = out_csv,
      stats = parc_df,
      note = msg,
      atlas_labels_path = paste0(atlas_labels_base, ".nii.gz"),
      roi_dir = roi_dir
    ))
  }
  
  
  
  # Define the reactive value to control the play state at the top level
  play_state <- reactiveVal(FALSE)
  
  # Dynamically detect mounted drives in /mnt/
  available_drives <- list.dirs("/mnt", recursive = FALSE, full.names = TRUE)
  
  # Extract the drive letters (last part of the path)
  drive_names <- basename(available_drives)
  
  # Create the volumes list dynamically with proper names
  volumes <- c(Home = "/home", Root = "/", setNames(available_drives, paste0("Windows_", toupper(drive_names))))
  
  
  # File chooser logic for PET Scan Input Path
  shinyFileChoose(input, 'btn_pet', roots = volumes, session = session)
  observeEvent(input$btn_pet, {
    if (!is.null(input$btn_pet)) {
      files <- shinyFiles::parseFilePaths(volumes, input$btn_pet)
      if (nrow(files) > 0) {
        selected_paths <- as.character(files$datapath)
        updateTextInput(session, "pet_input", value = paste(selected_paths, collapse = ","))
      }
    }
  })
  
  # File chooser logic for MRI Scans Input Path
  shinyFileChoose(input, 'btn_mri', roots = volumes, session = session)
  observeEvent(input$btn_mri, {
    if (!is.null(input$btn_mri)) {
      files <- shinyFiles::parseFilePaths(volumes, input$btn_mri)
      if (nrow(files) > 0) {
        selected_paths <- as.character(files$datapath)
        updateTextInput(session, "mri_input", value = paste(selected_paths, collapse = ","))
      }
    }
  })
  
  # File chooser logic for CSV Input Path
  shinyFileChoose(input, 'btn_csv', roots = volumes, session = session, filetypes = c('csv', 'xlsx', 'xls'))
  observeEvent(input$btn_csv, {
    if (!is.null(input$btn_csv)) {
      files <- shinyFiles::parseFilePaths(volumes, input$btn_csv)
      if (nrow(files) > 0) {
        selected_paths <- as.character(files$datapath)
        updateTextInput(session, "csv_input", value = selected_paths)
      }
    }
  })
  
  # File chooser logic for Control Region Mask Input Path
  shinyFileChoose(input, 'btn_control_mask', roots = volumes, session = session)
  observeEvent(input$btn_control_mask, {
    if (!is.null(input$btn_control_mask)) {
      files <- shinyFiles::parseFilePaths(volumes, input$btn_control_mask)
      if (nrow(files) > 0) {
        selected_paths <- as.character(files$datapath)
        updateTextInput(session, "control_mask_input", value = selected_paths)
      }
    }
  })
  
  # File chooser logic for Brain Region Mask Input Path
  shinyFileChoose(input, 'btn_brain_mask', roots = volumes, session = session)
  observeEvent(input$btn_brain_mask, {
    if (!is.null(input$btn_brain_mask)) {
      files <- shinyFiles::parseFilePaths(volumes, input$btn_brain_mask)
      if (nrow(files) > 0) {
        selected_paths <- as.character(files$datapath)
        updateTextInput(session, "brain_mask_input", value = selected_paths)
      }
    }
  })
  
  observeEvent(input$artifact_viewer, {
    # Path to the Dash app script
    dash_app_path <- file.path(getwd(), "Dashboard/ImagePropertiesDashboard.py")
    
    # Find an available port dynamically
    port <- as.integer(sample(49152:65535, 1))  # Use dynamic ports in the ephemeral range
    
    # Start the Dash app with the dynamic port and bind to 0.0.0.0
    # Remove 'nohup' for now to see any errors in the R console
    system(paste("python3", shQuote(dash_app_path), port, "0.0.0.0"), wait = FALSE)
    
    # Wait briefly to allow the Dash app to start
    Sys.sleep(2)
    
    # Construct the URL with the dynamic port
    url <- paste0("http://localhost:", port)
    
    # Open the Dash app in the Windows default browser
    system(paste("cmd.exe /C start", url))
  })
  
  
  # Create reactive values to store outputs
  reactive_values <- reactiveValues(output_choices = NULL)
  
  # Observe the preprocess button click
  observeEvent(input$preprocess, {
    
    withProgress(message = "Processing...", value = 0, {
      incProgress(0.02, detail = "Please wait")
      
      
      # Reset reactive values before any processing
      reactive_values$output_choices <- NULL
      
      # Reset the plot output and slider input
      output$brain_plot <- NULL
      updateSliderInput(session, "slice_number", min = 1, max = 1, value = 1)
      play_state(FALSE)
      updateActionButton(session, "play_slices", disabled = FALSE)
      updateActionButton(session, "stop_slices", disabled = TRUE)
      
      # Reset the plot and the dropdown choices
      output$brain_plot <- renderPlot(NULL)  # Ensure the plot is cleared
      updateSelectInput(session, "selected_output", choices = NULL)  # Clear the dropdown menu
      
      
      # Check if no preprocessing steps were selected
      if (length(input$preprocessing_steps) == 0) {
        showModal(modalDialog(
          title = "No Preprocessing Steps Selected",
          "Please select at least one preprocessing step to proceed.",
          easyClose = TRUE,
          footer = modalButton("OK")
        ))
        return()  # Stop further execution
      }
      
      # Function to remove file extensions
      remove_extensions <- function(file) {
        file <- gsub("\\.gz$", "", file)
        file <- gsub("\\.nii$", "", file)
        return(file)
      }
      
      # Process MRI Scan if brain extraction is selected
      if ("brain_extraction" %in% input$preprocessing_steps) {
        
        # Forcefully reset the plot and dropdown before starting brain extraction
        output$brain_plot <- renderPlot({})  # Forcefully clear the plot output
        updateSelectInput(session, "selected_output", choices = character(0))  # Clear dropdown menu completely
        
        
        # Check if MRI paths are provided
        if (nzchar(input$mri_input)) {
          mri_files <- unlist(strsplit(input$mri_input, ",")) # Splitting the input paths correctly
          
          # Initialize lists to store input and output file paths
          input_files <- list()
          output_files <- list()
          
          for (mri_file in mri_files) {
            # Normalize each MRI file path
            file_path <- normalizePath(mri_file, mustWork = TRUE)
            
            # Update the console messages to indicate which file is being processed
            consoleMessages(paste0(consoleMessages(), "Performing Brain Extraction on MRI scan: ", basename(file_path), "...\n"))
            
            # Define the output directory
            base_output_dir <- file.path(getwd(), "Output_Brain_Extraction")
            dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Create a subdirectory with the name of the file without extension
            file_name_without_extension <- remove_extensions(tools::file_path_sans_ext(basename(file_path)))
            output_dir <- file.path(base_output_dir, file_name_without_extension)
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Set output file path
            output_file <- file.path(output_dir, paste0(file_name_without_extension, ".nii.gz"))
            
            # Call the function for brain extraction using extrantsr
            extrantsr::fslbet_robust(file_path, output_file, remover = 'double_remove_neck')
            
            
            # Update console messages
            consoleMessages(paste0(consoleMessages(), "Brain extraction completed and saved to ", output_file, "\n"))
            
            # After brain extraction completes
            # Load images (for metrics + mid-slice PNGs inside report)
            input_img  <- readNIfTI(file_path, reorient = FALSE)
            output_img <- readNIfTI(output_file, reorient = FALSE)
            
            report_html <- generate_processing_report(
              step = "Brain Extraction",
              output_dir = output_dir,
              input_paths = list(MRI = file_path),
              output_paths = list(SkullStrippedMRI = output_file),
              parameters = list(method="extrantsr::fslbet_robust", remover="double_remove_neck"),
              logs = consoleMessages(),
              input_img = input_img,
              output_img = output_img
            )
            
            consoleMessages(paste0(consoleMessages(),
                                   "Report saved: ", report_html, "\n"
            ))
            
            # Store input and output file paths
            input_files[[basename(output_file)]] <- file_path
            output_files[[basename(output_file)]] <- output_file
          }
          
          # Update the selectInput with processed outputs
          updateSelectInput(session, "selected_output",
                            choices = names(output_files),
                            selected = names(output_files)[1])
          
          # Initialize reactive values to store the images
          mri_img <- reactiveVal()
          brain_img <- reactiveVal()
          
          # Load the images based on the selected output
          observeEvent(input$selected_output, {
            selected_name <- input$selected_output
            req(selected_name)
            
            # Load the MRI and brain-extracted images
            mri_file <- input_files[[selected_name]]
            brain_file <- output_files[[selected_name]]
            
            mri_img(readNIfTI(mri_file, reorient = FALSE))
            brain_img(readNIfTI(brain_file, reorient = FALSE))
            
            # Update slider input max value based on MRI dimensions
            img_dims <- dim(mri_img())
            updateSliderInput(session, "slice_number", min = 1, max = img_dims[3], value = round(img_dims[3] / 2))
          })
          
          # Load images for the first selection
          first_selection <- names(output_files)[1]
          if (!is.null(first_selection)) {
            selected_name <- first_selection
            mri_file <- input_files[[selected_name]]
            brain_file <- output_files[[selected_name]]
            
            mri_img(readNIfTI(mri_file, reorient = FALSE))
            brain_img(readNIfTI(brain_file, reorient = FALSE))
            
            # Update slider input max value based on MRI dimensions
            img_dims <- dim(mri_img())
            updateSliderInput(session, "slice_number", min = 1, max = img_dims[3], value = round(img_dims[3] / 2))
          }
          
          # Observe Previous Slice button click
          observeEvent(input$prev_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice > 1) {
              updateSliderInput(session, "slice_number", value = current_slice - 1)
            }
          })
          
          # Observe Next Slice button click
          observeEvent(input$next_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice < dim(mri_img())[3]) {
              updateSliderInput(session, "slice_number", value = current_slice + 1)
            }
          })
          
          # Observe Play button click
          observeEvent(input$play_slices, {
            play_state(TRUE)
            updateActionButton(session, "play_slices", disabled = TRUE)
            updateActionButton(session, "stop_slices", disabled = FALSE)
          })
          
          # Observe Stop button click
          observeEvent(input$stop_slices, {
            play_state(FALSE)
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
          })
          
          # Automatically update the slice number during play state
          observe({
            if (play_state()) {
              current_slice <- input$slice_number
              if (current_slice < dim(mri_img())[3]) {
                updateSliderInput(session, "slice_number", value = current_slice + 1)
                invalidateLater(100, session)  # Adjust the interval as needed
              } else {
                play_state(FALSE) # Stop play when the last slice is reached
                updateActionButton(session, "play_slices", disabled = FALSE)
                updateActionButton(session, "stop_slices", disabled = TRUE)
              }
            }
          })
          
          # Load and visualize the MRI and brain-extracted images
          output$brain_plot <- renderPlot({
            req(mri_img(), brain_img())
            slice_num <- as.integer(input$slice_number)
            
            # Set up the plotting area to have 1 row and 2 columns with a black background
            par(mfrow = c(1, 2), oma = c(2, 1, 2, 2), mar = c(2, 2, 2, 2), bg = "black", fg = "white")  # Adjust margins for titles
            
            # Function to plot a single slice with black background if empty
            plot_slice <- function(slice, title) {
              if (all(slice == 0)) {
                par(bg = "black")  # Set the background to black
                plot.new()
                plot.window(xlim = c(0, 1), ylim = c(0, 1))
                rect(0, 0, 1, 1, col = "black", border = "black")
                title(title, col.main = "white")
              } else {
                image(slice, col = gray(0:64/64), axes = FALSE, main = "", asp = 1)
                title(title, col.main = "white")
              }
            }
            
            # Display the selected slice for both images
            plot_slice(mri_img()[,,slice_num], "Original Input MRI")
            plot_slice(brain_img()[,,slice_num], "Skull-Stripped MRI")
            
            # Adjust position for the entire title block
            par(fig = c(0.5, 1, 0, 0.15), new = TRUE, bg = "black")
            
            # Add filename and titles at the bottom right corner of the plot
            #selected_name <- input$selected_output
            #mtext(paste("Input Filename:", selected_name), side = 1, line = -4, adj = 1, cex = 1.5, col = "white")
            
            
          }, res = 200)
          # Reset all inputs after brain extraction is complete
          updateTextInput(session, "pet_input", value = "")
          updateTextInput(session, "mri_input", value = "")
          updateTextInput(session, "csv_input", value = "")
          updateTextInput(session, "control_mask_input", value = "")
          updateTextInput(session, "brain_mask_input", value = "")
        } else {
          consoleMessages(paste0(consoleMessages(), "Please select MRI file(s) for Brain Extraction.\n"))
        }
      }
      
      
      # Process MRI Scan if MNI registration is selected
      if ("mni_registration" %in% input$preprocessing_steps) {
        
        # Forcefully reset the plot and dropdown before starting brain extraction
        output$brain_plot <- renderPlot({})  # Forcefully clear the plot output
        updateSelectInput(session, "selected_output", choices = character(0))  # Clear dropdown menu completely
        
        
        if (!is.null(input$mri_input) && nzchar(input$mri_input)) {  # Check if input$mri_input is not null and is a valid character
          mri_files <- unlist(strsplit(as.character(input$mri_input), ","))  # Splitting the input paths correctly
          consoleMessages(paste0(consoleMessages(), "Performing MNI Template Registration on multiple files...\n"))
          
          template_file <- file.path(getwd(), "templates", input$template_selection)
          
          # Initialize lists to store input and output file paths
          input_files <- list()
          output_files <- list()
          
          for (mri_file in mri_files) {
            file_path <- normalizePath(mri_file, mustWork = TRUE)
            
            # Define and create the output directory for each file
            base_output_dir <- file.path(getwd(), "Output_MNI_Template_Registration")
            file_name_without_extension <- remove_extensions(tools::file_path_sans_ext(basename(file_path)))
            mni_output_dir <- file.path(base_output_dir, file_name_without_extension)
            dir.create(mni_output_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Define paths for the brain-extracted MRI and the template
            mri_mask <- file.path(mni_output_dir, paste0(file_name_without_extension, "_brain.nii.gz"))
            template_mask <- file.path(mni_output_dir, "stripped_template.nii.gz")
            
            # Perform robust brain extraction using extrantsr
            extrantsr::fslbet_robust(file_path, mri_mask, remover = 'double_remove_neck')
            extrantsr::fslbet_robust(template_file, template_mask, remover = 'double_remove_neck')
            
            # Perform registration using FLIRT
            mri_registered <- file.path(mni_output_dir, paste0(file_name_without_extension, "_registered.nii.gz"))
            mri_mat_file <- file.path(mni_output_dir, "mri.mat")
            
            system2("flirt", args = c("-in", mri_mask,
                                      "-ref", template_mask,
                                      "-out", mri_registered,
                                      "-omat", mri_mat_file,
                                      "-dof", 12), stdout = TRUE)
            
            if (!file.exists(mri_registered)) {
              stop("Registration failed: output file not created")
            }
            
            consoleMessages(paste0(consoleMessages(), "MNI Template Registration completed for ", basename(file_path), " and saved to ", mni_output_dir, "\n"))
            
            # Load images for metrics
            reg_img <- readNIfTI(mri_registered, reorient = FALSE)
            ref_img <- readNIfTI(template_mask, reorient = FALSE)   # skull-stripped template
            in_img  <- readNIfTI(mri_mask, reorient = FALSE)        # skull-stripped input MRI
            
            
            # Optional parcellation (only if enabled)
            parc_res <- NULL
            if (isTRUE(input$do_parcellation)) {
              parc_res <- run_parcellation(
                target_img_path = mri_mask,
                output_dir = mni_output_dir,
                prefix = paste0(file_name_without_extension, "_mri_space"),
                mri_to_mni_mat = mri_mat_file,
                pet_to_mri_mat = NULL,
                pet_ref_path   = NULL,
                mri_ref_path   = mri_mask
              )
            }
            
            report_html <- generate_processing_report(
              step = "MNI Template Registration",
              output_dir = mni_output_dir,
              input_paths = list(SkullStrippedMRI = mri_mask, SkullStrippedTemplate = template_mask),
              output_paths = list(RegisteredMRI = mri_registered, AffineMat = mri_mat_file),
              parameters = list(dof=12, tool="flirt"),
              parcellation = parc_res,
              logs = consoleMessages(),
              input_img = in_img,
              output_img = reg_img,
              ref_img = ref_img,
              transform_mat_path = mri_mat_file
            )
            
            consoleMessages(paste0(consoleMessages(),
                                   "Report saved: ", report_html, "\n"
            ))
            
            
            
            # Store input and output file paths
            input_files[[basename(mri_registered)]] <- mri_mask
            output_files[[basename(mri_registered)]] <- mri_registered
          }
          
          # Update the selectInput with processed outputs
          updateSelectInput(session, "selected_output",
                            choices = names(output_files),
                            selected = names(output_files)[1])
          
          # Initialize reactive values to store the images
          mri_img <- reactiveVal()
          template_img <- reactiveVal()
          registered_img <- reactiveVal()
          
          # Load the images based on the selected output
          observeEvent(input$selected_output, {
            selected_name <- input$selected_output
            req(selected_name)
            
            # Load the MRI and registered images
            mri_file <- input_files[[selected_name]]
            registered_file <- output_files[[selected_name]]
            
            mri_img(readNIfTI(mri_file, reorient = FALSE))
            registered_img(readNIfTI(registered_file, reorient = FALSE))
            template_img(readNIfTI(template_mask, reorient = FALSE))
            
            # Calculate the minimum number of slices among the images
            max_slices <- min(dim(mri_img())[3], dim(template_img())[3], dim(registered_img())[3])
            
            # Update slider input max value based on minimum number of slices
            updateSliderInput(session, "slice_number", min = 1, max = max_slices, value = round(max_slices / 2))
          })
          
          # Load images for the first selection
          first_selection <- names(output_files)[1]
          if (!is.null(first_selection)) {
            selected_name <- first_selection
            mri_file <- input_files[[selected_name]]
            registered_file <- output_files[[selected_name]]
            
            mri_img(readNIfTI(mri_file, reorient = FALSE))
            registered_img(readNIfTI(registered_file, reorient = FALSE))
            template_img(readNIfTI(template_mask, reorient = FALSE))
            
            # Calculate the minimum number of slices among the images
            max_slices <- min(dim(mri_img())[3], dim(template_img())[3], dim(registered_img())[3])
            
            # Update slider input max value based on minimum number of slices
            updateSliderInput(session, "slice_number", min = 1, max = max_slices, value = round(max_slices / 2))
          }
          
          # Observe Previous Slice button click
          observeEvent(input$prev_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice > 1) {
              updateSliderInput(session, "slice_number", value = current_slice - 1)
            }
          })
          
          # Observe Next Slice button click
          observeEvent(input$next_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice < max_slices) {
              updateSliderInput(session, "slice_number", value = current_slice + 1)
            }
          })
          
          # Observe Play button click
          observeEvent(input$play_slices, {
            play_state(TRUE)
            updateActionButton(session, "play_slices", disabled = TRUE)
            updateActionButton(session, "stop_slices", disabled = FALSE)
          })
          
          # Observe Stop button click
          observeEvent(input$stop_slices, {
            play_state(FALSE)
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
          })
          
          # Automatically update the slice number during play state
          observe({
            if (play_state()) {
              current_slice <- input$slice_number
              if (current_slice < max_slices) {
                updateSliderInput(session, "slice_number", value = current_slice + 1)
                invalidateLater(100, session)  # Adjust the interval as needed
              } else {
                play_state(FALSE) # Stop play when the last slice is reached
                updateActionButton(session, "play_slices", disabled = FALSE)
                updateActionButton(session, "stop_slices", disabled = TRUE)
              }
            }
          })
          
          # Render the plot
          output$brain_plot <- renderPlot({
            req(mri_img(), template_img(), registered_img())
            slice_num <- as.integer(input$slice_number)
            
            # Set up the plotting area to have 1 row and 3 columns with a black background
            par(mfrow = c(1, 3), oma = c(2, 1, 2, 2), mar = c(2, 2, 2, 2), bg = "black", fg = "white")
            
            # Function to plot a single slice with black background if empty
            plot_slice <- function(img, slice_num, title) {
              if (slice_num <= dim(img)[3]) {
                slice <- img[,,slice_num]
                if (all(slice == 0)) {
                  par(bg = "black")
                  plot.new()
                  plot.window(xlim = c(0, 1), ylim = c(0, 1))
                  rect(0, 0, 1, 1, col = "black", border = "black")
                  title(title, col.main = "white")
                } else {
                  image(slice, col = gray(0:64/64), axes = FALSE, main = "", asp = 1)
                  title(title, col.main = "white")
                }
              } else {
                par(bg = "black")
                plot.new()
                plot.window(xlim = c(0, 1), ylim = c(0, 1))
                rect(0, 0, 1, 1, col = "black", border = "black")
                title(title, col.main = "white")
              }
            }
            
            # Display the selected slice for all three images
            plot_slice(mri_img(), slice_num, "Skull-Stripped Input MRI\n(native space)")
            plot_slice(template_img(), slice_num, "Skull-Stripped MNI Template\n(output space)")
            plot_slice(registered_img(), slice_num, "Registered MRI\n(output space)")
            
            # Add filename at the bottom
            #selected_name <- input$selected_output
            #mtext(paste("Input Filename:", selected_name), side = 1, line = -2, adj = 0.5, cex = 1.5, col = "white", outer = TRUE)
            
            
          }, res = 200)
          # Reset all inputs after brain extraction is complete
          updateTextInput(session, "pet_input", value = "")
          updateTextInput(session, "mri_input", value = "")
          updateTextInput(session, "csv_input", value = "")
          updateTextInput(session, "control_mask_input", value = "")
          updateTextInput(session, "brain_mask_input", value = "")
          
        } else {
          consoleMessages(paste0(consoleMessages(), "Please select MRI file(s) for MNI Template Registration.\n"))
        }
      }
      
      
      
      # Process MRI Scan if MRI-PET Co-registration is selected
      if ("mri_coregistration" %in% input$preprocessing_steps) {
        
        # Forcefully reset the plot and dropdown before starting brain extraction
        output$brain_plot <- renderPlot({})  # Forcefully clear the plot output
        updateSelectInput(session, "selected_output", choices = character(0))  # Clear dropdown menu completely
        
        
        if (!is.null(input$csv_input) && nzchar(input$csv_input)) {  # Added check to ensure csv_input is not NULL and is a non-empty string
          consoleMessages(paste0(consoleMessages(), "Performing MRI-PET Co-registration using paths from CSV/Excel...\n"))
          
          # Read the CSV or Excel file
          csv_file <- as.character(input$csv_input)
          
          # Check file extension
          file_extension <- tools::file_ext(csv_file)
          
          # Read data based on file extension
          if (file_extension == "csv") {
            df <- tryCatch({
              read.csv(csv_file, stringsAsFactors = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading CSV file: ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(df)) return()
          } else if (file_extension %in% c("xls", "xlsx")) {
            df <- tryCatch({
              read_excel(csv_file)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading Excel file: ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(df)) return()
          } else {
            consoleMessages(paste0(consoleMessages(), "Unsupported file type. Please upload a CSV or Excel file.\n"))
            return()
          }
          
          # Check if required columns are present
          if (!all(c('MRI_InputPath', 'PET_InputPath') %in% names(df))) {
            consoleMessages(paste0(consoleMessages(), "CSV/Excel file must contain columns 'MRI_InputPath' and 'PET_InputPath'.\n"))
            return()
          }
          
          # Normalize path separators
          df$MRI_InputPath <- gsub("\\\\", "/", df$MRI_InputPath)
          df$PET_InputPath <- gsub("\\\\", "/", df$PET_InputPath)
          
          # Initialize lists to store input and output file paths
          input_files <- list()
          output_files <- list()
          
          # Prepare the template
          template_file <- file.path(getwd(), "templates", input$template_selection)
          temp.img <- tryCatch({
            readNIfTI(template_file, reorient = FALSE)
          }, error = function(e) {
            consoleMessages(paste0(consoleMessages(), "Error reading template file: ", e$message, "\n"))
            return(NULL)
          })
          if (is.null(temp.img)) return()
          
          ss_template <- tryCatch({
            extrantsr::fslbet_robust(temp.img, remover = 'double_remove_neck')
          }, error = function(e) {
            consoleMessages(paste0(consoleMessages(), "Error performing brain extraction on template: ", e$message, "\n"))
            return(NULL)
          })
          if (is.null(ss_template)) return()
          
          # Process each MRI-PET pair from the CSV/Excel file
          for (i in 1:nrow(df)) {
            mri_file <- as.character(df$MRI_InputPath[i])
            pet_file <- as.character(df$PET_InputPath[i])
            
            # Normalize path separators again (in case)
            mri_file <- gsub("\\\\", "/", mri_file)
            pet_file <- gsub("\\\\", "/", pet_file)
            
            # Check if the files exist
            if (!file.exists(mri_file) || !file.exists(pet_file)) {
              consoleMessages(paste0(consoleMessages(), "File not found: ", mri_file, " or ", pet_file, "\n"))
              next  # Skip to the next pair
            }
            
            # Read MRI and PET images with reorient = FALSE
            mri.img <- tryCatch({
              readNIfTI(mri_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading MRI file (Pair ", i, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(mri.img)) { next }
            
            pet.img <- tryCatch({
              readNIfTI(pet_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading PET file (Pair ", i, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(pet.img)) { next }
            
            # Perform brain extraction on MRI
            mri.mask <- tryCatch({
              extrantsr::fslbet_robust(mri.img, remover = 'double_remove_neck')
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error performing brain extraction on MRI (Pair ", i, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(mri.mask)) { next }
            
            # Define the output directory
            base_output_dir <- file.path(getwd(), "Output_MRI_PET_Co-Registration")
            dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Create a subdirectory with the name of the PET file without extension
            file_name_without_extension <- tryCatch({
              remove_extensions(tools::file_path_sans_ext(basename(pet_file)))
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error processing PET file name (Pair ", i, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(file_name_without_extension)) { next }
            
            mni_output_dir <- file.path(base_output_dir, file_name_without_extension)
            dir.create(mni_output_dir, recursive = TRUE, showWarnings = FALSE)
            
            # Register the original PET to the original MRI using FLIRT
            pet2mri_mat_file <- file.path(mni_output_dir, "pet2mri.mat")
            registered_pet <- tryCatch({
              flirt(infile = pet.img, reffile = mri.img, dof = 12, omat = pet2mri_mat_file)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error during PET to MRI registration (Pair ", i, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(registered_pet)) { next }
            
            # Register the skull-stripped MRI to the skull-stripped template using FLIRT
            mri_mat_file <- file.path(mni_output_dir, "mri.mat")
            mri_registered <- tryCatch({
              flirt(infile = mri.mask, reffile = ss_template, omat = mri_mat_file, dof = 12)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error during MRI to Template registration (Pair ", i, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(mri_registered)) { next }
            
            # Apply the transformation from MRI registration to the PET image using FLIRT
            transformed_pet_file <- file.path(mni_output_dir, paste0("MRI-PET_Coregistered_", basename(pet_file)))
            #transformed_pet_file <- file.path(mni_output_dir, "transformed_pet.nii.gz")
            
            # Apply resulting mask to MRI registered PET (binary mask)
            mask <- mri.mask > 0
            pet.mask <- tryCatch({
              fslmask(registered_pet, mask = mask)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error applying mask (Pair ", i, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(pet.mask)) { next }
            
            # Apply transformation with error handling
            flirt_apply_result <- tryCatch({
              flirt_apply(infile = pet.mask, reffile = ss_template, initmat = mri_mat_file, outfile = transformed_pet_file)
            }, error = function(e) {
              if (grepl("Transformation is not simple, cannot reorient!", e$message, ignore.case = TRUE)) {
                consoleMessages(paste0(consoleMessages(), "Transformation is not simple for pair ", i, ". Skipping reorientation.\n"))
              } else {
                consoleMessages(paste0(consoleMessages(), "Error applying transformation (Pair ", i, "): ", e$message, "\n"))
              }
              return(NULL)
            })
            if (is.null(flirt_apply_result)) { next }
            
            consoleMessages(paste0(consoleMessages(), "MRI-PET Co-registration completed for pair ", i, " and saved to ", transformed_pet_file, "\n"))
            
            # Load images for metrics
            mri_viz <- readNIfTI(mri_file, reorient = FALSE)
            pet_out <- readNIfTI(transformed_pet_file, reorient = FALSE)
            
            # Optional parcellation (only if enabled)
            parc_res <- NULL
            if (isTRUE(input$do_parcellation) && ("mri_coregistration" %in% input$preprocessing_steps)) {
              parc_res <- run_parcellation(
                target_img_path = transformed_pet_file,   # keep your existing target (doesn't break report)
                output_dir = mni_output_dir,
                prefix = paste0("pair_", i, "_pet_space"),
                mri_to_mni_mat = mri_mat_file,            # MRI->MNI from your existing template registration
                pet_to_mri_mat = pet2mri_mat_file,        # PET->MRI we just saved
                pet_ref_path   = pet_file,                # ORIGINAL PET defines PET space
                mri_ref_path   = mri_file                 # fallback only
              )
            }
            
            
            report_html <- generate_processing_report(
              step = "MRI–PET Co-registration",
              output_dir = mni_output_dir,
              input_paths = list(MRI = mri_file, PET = pet_file),
              output_paths = list(CoregPET = transformed_pet_file, MRIAffine = mri_mat_file),
              parameters = list(dof=12, tool="flirt / fslr::flirt_apply"),
              parcellation = parc_res,
              logs = consoleMessages(),
              input_img = mri_viz,
              output_img = pet_out,
              ref_img = ss_template,              # optional: you’re outputting in template space
              transform_mat_path = mri_mat_file
            )
            
            consoleMessages(paste0(consoleMessages(),
                                   "Report saved: ", report_html, "\n"
            ))
            
            
            # Store input and output file paths
            input_files[[basename(transformed_pet_file)]] <- pet_file
            output_files[[basename(transformed_pet_file)]] <- transformed_pet_file
          }
          
          # Update the selectInput with processed outputs
          updateSelectInput(session, "selected_output",
                            choices = names(output_files),
                            selected = names(output_files)[1])
          
          # Initialize reactive values to store the images
          pet_img <- reactiveVal()
          ref_img_viz <- reactiveVal()
          pet_reg_viz <- reactiveVal()
          
          # Load the images based on the selected output
          observeEvent(input$selected_output, {
            selected_name <- input$selected_output
            req(selected_name)
            
            # Load the original PET and co-registered PET images
            pet_file <- input_files[[selected_name]]
            transformed_pet_file <- output_files[[selected_name]]
            
            pet_img(readNIfTI(pet_file, reorient = FALSE))
            ref_img_viz(ss_template)
            pet_reg_viz(readNIfTI(transformed_pet_file, reorient = FALSE))
            
            # Slider should follow the output/reference space
            max_slices <- min(dim(ref_img_viz())[3], dim(pet_reg_viz())[3])
            updateSliderInput(session, "slice_number", min = 1, max = max_slices, value = round(max_slices / 2))
          })
          
          # Load images for the first selection
          first_selection <- names(output_files)[1]
          if (!is.null(first_selection)) {
            selected_name <- first_selection
            pet_file <- input_files[[selected_name]]
            transformed_pet_file <- output_files[[selected_name]]
            
            pet_img(readNIfTI(pet_file, reorient = FALSE))
            ref_img_viz(ss_template)
            pet_reg_viz(readNIfTI(transformed_pet_file, reorient = FALSE))
            
            # Slider should follow the output/reference space
            max_slices <- min(dim(ref_img_viz())[3], dim(pet_reg_viz())[3])
            updateSliderInput(session, "slice_number", min = 1, max = max_slices, value = round(max_slices / 2))
          }
          
          # Observe Previous Slice button click
          observeEvent(input$prev_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice > 1) {
              updateSliderInput(session, "slice_number", value = current_slice - 1)
            }
          })
          
          # Observe Next Slice button click
          observeEvent(input$next_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice < dim(pet_reg_viz())[3]) {
              updateSliderInput(session, "slice_number", value = current_slice + 1)
            }
          })
          
          # Observe Play button click
          observeEvent(input$play_slices, {
            play_state(TRUE)
            updateActionButton(session, "play_slices", disabled = TRUE)
            updateActionButton(session, "stop_slices", disabled = FALSE)
          })
          
          # Observe Stop button click
          observeEvent(input$stop_slices, {
            play_state(FALSE)
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
          })
          
          # Automatically update the slice number during play state
          observe({
            if (play_state()) {
              current_slice <- input$slice_number
              if (current_slice < dim(pet_reg_viz())[3]) {
                updateSliderInput(session, "slice_number", value = current_slice + 1)
                invalidateLater(100, session)  # Adjust the interval as needed
              } else {
                play_state(FALSE) # Stop play when the last slice is reached
                updateActionButton(session, "play_slices", disabled = FALSE)
                updateActionButton(session, "stop_slices", disabled = TRUE)
              }
            }
          })
          
          # Render the plot
          # Render the plot
          output$brain_plot <- renderPlot({
            req(pet_img(), ref_img_viz(), pet_reg_viz())
            
            # Slice in output/reference space
            ref_slice_num <- as.integer(input$slice_number)
            
            # Map reference-space slice position to original PET slice position
            ref_nz  <- dim(ref_img_viz())[3]
            orig_nz <- dim(pet_img())[3]
            
            if (ref_nz > 1 && orig_nz > 1) {
              rel_pos <- (ref_slice_num - 1) / (ref_nz - 1)
              orig_slice_num <- round(rel_pos * (orig_nz - 1)) + 1
            } else {
              orig_slice_num <- 1
            }
            
            # Safety clamp
            orig_slice_num <- max(1, min(orig_slice_num, orig_nz))
            ref_slice_num  <- max(1, min(ref_slice_num, ref_nz))
            
            # Set up the plotting area to have 1 row and 3 columns with a black background
            par(mfrow = c(1, 3), oma = c(2, 1, 2, 2), mar = c(2, 2, 2, 2), bg = "black", fg = "white")
            
            # Function to plot a single slice with black background if empty
            plot_slice <- function(img, slice_num, title) {
              if (slice_num <= dim(img)[3]) {
                slice <- img[,,slice_num]
                if (all(slice == 0)) {
                  par(bg = "black")
                  plot.new()
                  plot.window(xlim = c(0, 1), ylim = c(0, 1))
                  rect(0, 0, 1, 1, col = "black", border = "black")
                  title(title, col.main = "white")
                } else {
                  image(slice, col = gray(0:64/64), axes = FALSE, main = "", asp = 1)
                  title(title, col.main = "white")
                }
              } else {
                par(bg = "black")
                plot.new()
                plot.window(xlim = c(0, 1), ylim = c(0, 1))
                rect(0, 0, 1, 1, col = "black", border = "black")
                title(title, col.main = "white")
              }
            }
            
            # Left: original PET in native space, using proportionally matched slice
            plot_slice(
              pet_img(),
              orig_slice_num,
              paste0("Original PET Image\n(native space)")
            )
            
            # Middle: reference image defining output space
            plot_slice(
              ref_img_viz(),
              ref_slice_num,
              paste0("Reference Template\n(output space)")
            )
            
            # Right: transformed PET in same space as middle panel
            plot_slice(
              pet_reg_viz(),
              ref_slice_num,
              paste0("Co-registered PET Image\n(output space)")
            )
            
          }, res = 200)
          # Reset all inputs after brain extraction is complete
          updateTextInput(session, "pet_input", value = "")
          updateTextInput(session, "mri_input", value = "")
          updateTextInput(session, "csv_input", value = "")
          updateTextInput(session, "control_mask_input", value = "")
          updateTextInput(session, "brain_mask_input", value = "")
          
        } else {
          consoleMessages(paste0(consoleMessages(), "Please provide a CSV/Excel file with MRI and PET file paths for Co-registration.\n"))
        }
      }
      
      
      
      # COMBAT Harmonization
      if ("combat_harmonization" %in% input$preprocessing_steps) {
        
        # Forcefully reset the plot and dropdown before starting brain extraction
        output$brain_plot <- renderPlot({})  # Forcefully clear the plot output
        updateSelectInput(session, "selected_output", choices = character(0))  # Clear dropdown menu completely
        
        
        # Ensure input$csv_input is not NULL and is a valid non-empty character string
        if (!is.null(input$csv_input) && nzchar(input$csv_input)) {
          csv_input <- as.character(input$csv_input)  # Ensure it is a character
          
          consoleMessages(paste0(consoleMessages(), "Performing COMBAT Harmonization...\n"))
          
          # Read covariate data
          covariate_file <- csv_input
          file_extension <- tools::file_ext(covariate_file)
          
          # Read data based on file extension with error handling
          covariates <- tryCatch({
            if (file_extension == "csv") {
              read.csv(covariate_file, stringsAsFactors = FALSE)
            } else if (file_extension %in% c("xls", "xlsx")) {
              read_excel(covariate_file)
            } else {
              stop("Unsupported covariate file type. Please upload a CSV or Excel file.")
            }
          }, error = function(e) {
            consoleMessages(paste0(consoleMessages(), "Error reading covariate file: ", e$message, "\n"))
            return(NULL)
          })
          if (is.null(covariates)) return()
          
          # Check if required columns are present
          if (!all(c('Filename', 'Batch') %in% names(covariates))) {
            consoleMessages(paste0(consoleMessages(), "Covariate CSV/Excel file must contain columns 'Filename' and 'Batch'.\n"))
            return()
          }
          
          # Normalize paths in covariates$Filename
          covariates$Filename <- gsub("\\\\", "/", covariates$Filename)
          covariates$Filename_normalized <- tryCatch({
            normalizePath(covariates$Filename, winslash = "/", mustWork = FALSE)
          }, error = function(e) {
            consoleMessages(paste0("Error normalizing file paths: ", e$message, "\n"))
            return(NULL)
          })
          
          # Identify missing covariate files
          missing_covariate_files <- covariates$Filename[!file.exists(covariates$Filename_normalized)]
          if (length(missing_covariate_files) > 0) {
            consoleMessages(paste0(
              consoleMessages(),
              "The following covariate files do not exist:\n",
              paste(missing_covariate_files, collapse = "\n"),
              "\nPlease check the file paths in the CSV/Excel file and try again.\n"
            ))
            return()
          }
          
          # Proceed with normalized covariate paths
          covariates$Filename <- covariates$Filename_normalized
          
          # Extract PET file paths from covariates
          pet_files <- covariates$Filename
          
          # Ensure pet_files is a character vector
          pet_files <- as.character(pet_files)
          
          # Ensure covariate order matches the image data
          image_filenames <- data.frame(Filename = pet_files, stringsAsFactors = FALSE)
          covariates <- merge(image_filenames, covariates, by = "Filename", all.x = TRUE)
          
          # Check for missing batch information
          if (any(is.na(covariates$Batch))) {
            consoleMessages(paste0(consoleMessages(), "Missing batch information for some files.\n"))
            return()
          }
          
          # Initialize lists to store image data and dimensions
          image_data_list <- list()
          image_dims <- list()
          
          # Initialize lists to store input and output file paths
          input_files <- list()
          output_files <- list()
          
          # Loop through PET files and read data with error handling
          for (i in seq_along(pet_files)) {
            pet_file <- pet_files[i]
            
            # Check if the PET file exists
            if (!file.exists(pet_file)) {
              consoleMessages(paste0(consoleMessages(), "PET file does not exist: ", pet_file, "\n"))
              next  # Skip to the next file
            }
            
            # Read PET image with error handling
            pet_img <- tryCatch({
              readNIfTI(pet_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading PET file (", pet_file, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(pet_img)) { next }
            
            image_dims[[i]] <- dim(pet_img)
            image_data_list[[i]] <- as.numeric(pet_img)
            
            # Store input file paths with basename including extension
            input_files[[basename(pet_file)]] <- pet_file
          }
          
          # Check if any images were successfully read
          if (length(image_data_list) == 0) {
            consoleMessages(paste0(consoleMessages(), "No PET images were successfully read. COMBAT Harmonization aborted.\n"))
            return()
          }
          
          # Convert list to matrix (features x subjects)
          image_data_matrix <- do.call(cbind, image_data_list)
          
          # Get batch information
          batch <- covariates$Batch
          
          # Check if batch and data align properly
          if (length(batch) != ncol(image_data_matrix)) {
            consoleMessages(paste0(consoleMessages(), "Mismatch between batch info and number of PET scans.\n"))
            return()
          }
          
          # Identify and remove constant features
          feature_variance <- apply(image_data_matrix, 1, var)
          non_constant_features <- feature_variance != 0
          constant_features <- feature_variance == 0
          
          if (all(constant_features)) {
            consoleMessages(paste0(consoleMessages(), "All features are constant across samples. Cannot perform COMBAT Harmonization.\n"))
            return()
          }
          
          # Filter out constant features
          image_data_matrix_filtered <- image_data_matrix[non_constant_features, ]
          
          # Perform COMBAT Harmonization on non-constant features with error handling
          combat_results <- tryCatch({
            neuroCombat(dat = image_data_matrix_filtered, batch = batch)
          }, error = function(e) {
            consoleMessages(paste0(consoleMessages(), "Error during COMBAT Harmonization: ", e$message, "\n"))
            return(NULL)
          })
          if (is.null(combat_results)) return()
          
          harmonized_data <- combat_results$dat.combat
          
          # Reconstruct the full data matrix including constant features
          harmonized_data_full <- matrix(0, nrow = nrow(image_data_matrix), ncol = ncol(image_data_matrix))
          harmonized_data_full[non_constant_features, ] <- harmonized_data
          harmonized_data_full[constant_features, ] <- image_data_matrix[constant_features, ]  # Retain original values for constant features
          
          # Save harmonized images with error handling
          base_output_dir <- file.path(getwd(), "Output_COMBAT_Harmonization")
          dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
          
          for (i in seq_along(pet_files)) {
            pet_file <- pet_files[i]
            harmonized_vector <- harmonized_data_full[, i]
            
            # Reshape the harmonized data back into the original 3D array structure
            harmonized_array <- tryCatch({
              array(harmonized_vector, dim = image_dims[[i]])
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reshaping harmonized data for file ", pet_file, ": ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(harmonized_array)) { next }
            
            # Read the original NIfTI file again to copy its header with error handling
            original_nifti <- tryCatch({
              readNIfTI(pet_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading original PET file for header information (", pet_file, "): ", e$message, "\n"))
              return(NULL)
            })
            if (is.null(original_nifti)) { next }
            
            # Create a new NIfTI object with the harmonized data
            harmonized_nifti <- original_nifti
            harmonized_nifti[] <- harmonized_array  # Replace the image data
            
            # Define the output file path without extension
            output_file <- file.path(
              base_output_dir,
              paste0(
                "Combat_harmonized_",                            
                sub("\\.nii\\.gz$", "", basename(pet_file))  # Remove .nii.gz
              )
            )
            
            # Save the harmonized NIfTI object as a file with error handling
            tryCatch({
              writeNIfTI(harmonized_nifti, filename = output_file)  # writeNIfTI appends .nii.gz
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error writing harmonized PET file (", output_file, "): ", e$message, "\n"))
              return(NULL)
            })
            
            # Store output file paths with the full filename including extension
            output_files[[paste0("Combat_harmonized_", basename(pet_file))]] <- paste0(output_file, ".nii.gz")
            
            # Debugging: Print output file path
            consoleMessages(paste0("Saved harmonized file: ", paste0(output_file, ".nii.gz"), "\n"))
          }
          
          consoleMessages(paste0(consoleMessages(), "COMBAT Harmonization completed. Outputs saved to ", base_output_dir, "\n"))
          
          # --- COMBAT report (per output file) ---
          in_img  <- readNIfTI(pet_file, reorient = FALSE)
          out_img <- readNIfTI(paste0(output_file, ".nii.gz"), reorient = FALSE)
          
          report_html <- generate_processing_report(
            step = "COMBAT Harmonization",
            output_dir = base_output_dir,   # OR better: create per-file folder, see note below
            input_paths = list(PET = pet_file, Covariates = covariate_file),
            output_paths = list(HarmonizedPET = paste0(output_file, ".nii.gz")),
            parameters = list(
              method = "neuroCombat",
              batch_column = "Batch",
              constant_features_retained = TRUE
            ),
            logs = consoleMessages(),
            input_img = in_img,
            output_img = out_img
          )
          
          consoleMessages(paste0(consoleMessages(), "Report saved: ", report_html, "\n"))
          
          # Update the selectInput with processed outputs
          updateSelectInput(session, "selected_output",
                            choices = names(output_files),
                            selected = names(output_files)[1])
          
          # Initialize reactive values to store the images
          input_image <- reactiveVal()
          harmonized_image <- reactiveVal()
          
          # Load the images based on the selected output
          observeEvent(input$selected_output, {
            selected_name <- as.character(input$selected_output)
            req(selected_name)
            
            # Debugging: Print selected output name
            consoleMessages(paste0("Selected output: ", selected_name, "\n"))
            
            # Correctly map back to the original pet_file by removing the prefix
            pet_file_key <- sub("^Combat_harmonized_", "", selected_name)
            
            # Debugging: Print mapped pet_file_key
            consoleMessages(paste0("Mapped pet_file_key: ", pet_file_key, "\n"))
            
            pet_file <- input_files[[pet_file_key]]
            output_file <- output_files[[selected_name]]
            
            # Debugging: Print pet_file and output_file
            consoleMessages(paste0("pet_file: ", pet_file, "\n"))
            consoleMessages(paste0("output_file: ", output_file, "\n"))
            
            # Check that the pet_file and output_file are not NULL
            if (is.null(pet_file) || is.null(output_file)) {
              consoleMessages(paste0(consoleMessages(), "Error: File paths for selected output are invalid or missing.\n"))
              return()
            }
            
            # Load NIfTI images with error handling
            input_img <- tryCatch({
              readNIfTI(pet_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading input PET file (", pet_file, "): ", e$message, "\n"))
              return(NULL)
            })
            harmonized_img <- tryCatch({
              readNIfTI(output_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading harmonized PET file (", output_file, "): ", e$message, "\n"))
              return(NULL)
            })
            
            # Assign images to reactive values if successfully read
            if (!is.null(input_img) && !is.null(harmonized_img)) {
              input_image(input_img)
              harmonized_image(harmonized_img)
              
              # Update slider input max value based on image dimensions
              img_dims <- dim(harmonized_image())
              updateSliderInput(session, "slice_number", min = 1, max = img_dims[3], value = round(img_dims[3] / 2))
            } else {
              consoleMessages(paste0(consoleMessages(), "Error: Failed to load images for plotting.\n"))
            }
          })
          
          # Load images for the first selection
          first_selection <- names(output_files)[1]
          if (!is.null(first_selection)) {
            selected_name <- first_selection
            pet_file_key <- sub("^Combat_harmonized_", "", selected_name)
            pet_file <- input_files[[pet_file_key]]
            output_file <- output_files[[selected_name]]
            
            # Debugging: Print initial mapping
            consoleMessages(paste0("Initial selection: ", selected_name, "\n"))
            consoleMessages(paste0("Initial pet_file_key: ", pet_file_key, "\n"))
            consoleMessages(paste0("Initial pet_file: ", pet_file, "\n"))
            consoleMessages(paste0("Initial output_file: ", output_file, "\n"))
            
            # Load NIfTI images with error handling
            input_img <- tryCatch({
              readNIfTI(pet_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading input PET file (", pet_file, "): ", e$message, "\n"))
              return(NULL)
            })
            harmonized_img <- tryCatch({
              readNIfTI(output_file, reorient = FALSE)
            }, error = function(e) {
              consoleMessages(paste0(consoleMessages(), "Error reading harmonized PET file (", output_file, "): ", e$message, "\n"))
              return(NULL)
            })
            
            # Assign images to reactive values if successfully read
            if (!is.null(input_img) && !is.null(harmonized_img)) {
              input_image(input_img)
              harmonized_image(harmonized_img)
              
              # Update slider input max value based on image dimensions
              img_dims <- dim(harmonized_image())
              updateSliderInput(session, "slice_number", min = 1, max = img_dims[3], value = round(img_dims[3] / 2))
            } else {
              consoleMessages(paste0(consoleMessages(), "Error: Failed to load images for initial plotting.\n"))
            }
          }
          
          # Observe Previous Slice button click
          observeEvent(input$prev_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice > 1) {
              updateSliderInput(session, "slice_number", value = current_slice - 1)
            }
          })
          
          # Observe Next Slice button click
          observeEvent(input$next_slice, {
            play_state(FALSE) # Stop play when manually changing slices
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
            current_slice <- input$slice_number
            if (current_slice < dim(harmonized_image())[3]) {
              updateSliderInput(session, "slice_number", value = current_slice + 1)
            }
          })
          
          # Observe Play button click
          observeEvent(input$play_slices, {
            play_state(TRUE)
            updateActionButton(session, "play_slices", disabled = TRUE)
            updateActionButton(session, "stop_slices", disabled = FALSE)
          })
          
          # Observe Stop button click
          observeEvent(input$stop_slices, {
            play_state(FALSE)
            updateActionButton(session, "play_slices", disabled = FALSE)
            updateActionButton(session, "stop_slices", disabled = TRUE)
          })
          
          # Initialize play_state reactiveVal if not already
          if (!exists("play_state")) {
            play_state <- reactiveVal(FALSE)
          }
          
          # Automatically update the slice number during play state
          observe({
            if (play_state()) {
              current_slice <- input$slice_number
              if (current_slice < dim(harmonized_image())[3]) {
                updateSliderInput(session, "slice_number", value = current_slice + 1)
                invalidateLater(100, session)  # Adjust the interval as needed
              } else {
                play_state(FALSE) # Stop play when the last slice is reached
                updateActionButton(session, "play_slices", disabled = FALSE)
                updateActionButton(session, "stop_slices", disabled = TRUE)
              }
            }
          })
          
          # Render the plot
          output$brain_plot <- renderPlot({
            req(input_image(), harmonized_image())  # Ensure both images are available
            slice_num <- as.integer(input$slice_number)
            
            # Ensure slice number is within valid range
            img_dims <- dim(harmonized_image())
            if (is.na(slice_num) || slice_num < 1 || slice_num > img_dims[3]) {
              consoleMessages(paste0(consoleMessages(), "Invalid slice number: ", slice_num, "\n"))
              return()
            }
            
            # Set up the plotting area with black background
            par(mfrow = c(1, 2), oma = c(2, 1, 2, 2), mar = c(2, 2, 2, 2), bg = "black", fg = "white")
            
            # Define a function to plot a single slice with black background if empty
            plot_slice <- function(img, slice_num, title_text) {
              if (slice_num <= dim(img)[3]) {
                slice <- img[,,slice_num]
                if (all(slice == 0)) {
                  # Handle empty slices with black background
                  plot.new()
                  plot.window(xlim = c(0, 1), ylim = c(0, 1))
                  rect(0, 0, 1, 1, col = "black", border = "black")
                  title(title_text, col.main = "white")
                } else {
                  # Plot the actual slice
                  image(t(apply(slice, 2, rev)), col = gray(0:64/64), axes = FALSE, main = "", asp = 1)
                  title(title_text, col.main = "white")
                }
              } else {
                # If slice number exceeds dimensions, plot black background
                plot.new()
                plot.window(xlim = c(0, 1), ylim = c(0, 1))
                rect(0, 0, 1, 1, col = "black", border = "black")
                title(title_text, col.main = "white")
              }
            }
            
            # Plot the selected slice for both images
            plot_slice(input_image(), slice_num, "Original PET Image")
            plot_slice(harmonized_image(), slice_num, "COMBAT Harmonized Image")
            
            # Add the filename at the bottom of the plot
            #selected_name <- input$selected_output
            #mtext(paste("Input Filename:", selected_name), side = 1, line = -2, adj = 0.5, cex = 1.5, col = "white", outer = TRUE)
            
          }, res = 200)
          # Reset all inputs after brain extraction is complete
          updateTextInput(session, "pet_input", value = "")
          updateTextInput(session, "mri_input", value = "")
          updateTextInput(session, "csv_input", value = "")
          updateTextInput(session, "control_mask_input", value = "")
          updateTextInput(session, "brain_mask_input", value = "")
          
          
        } else {
          consoleMessages(paste0(consoleMessages(), "Please provide a CSV/Excel file with PET image file path and their batch details for COMBAT Harmonization.\n"))
        }
      }
      
      
      
      
      
      # RAVEL Normalization
      if ("ravel_normalization" %in% input$preprocessing_steps) {
        
        # Forcefully reset the plot and dropdown before starting brain extraction
        output$brain_plot <- renderPlot({})  # Forcefully clear the plot output
        updateSelectInput(session, "selected_output", choices = character(0))  # Clear dropdown menu completely
        
        
        if (nzchar(input$pet_input) && nzchar(input$control_mask_input) &&
            nzchar(input$brain_mask_input)) {
          
          consoleMessages(paste0(consoleMessages(), "Performing RAVEL Normalization...\n"))
          
          # Read PET scan paths
          pet_files <- unlist(strsplit(input$pet_input, ","))
          pet_files <- normalizePath(pet_files, winslash = "/", mustWork = TRUE)
          
          # Prepare the masks
          control_mask_file <- normalizePath(input$control_mask_input, winslash = "/", mustWork = TRUE)
          brain_mask_file <- normalizePath(input$brain_mask_input, winslash = "/", mustWork = TRUE)
          
          # Read the brain and control masks as NIfTI objects
          brain_mask_nifti <- readNIfTI(brain_mask_file, reorient = FALSE)
          control_mask_nifti <- readNIfTI(control_mask_file, reorient = FALSE)
          
          # Ensure masks are binary
          brain_mask_nifti <- brain_mask_nifti > 0
          control_mask_nifti <- control_mask_nifti > 0
          
          # Check that the dimensions of the masks match the PET images
          first_pet_nifti <- readNIfTI(pet_files[1], reorient = FALSE)
          if (!all(dim(brain_mask_nifti) == dim(first_pet_nifti))) {
            consoleMessages("Error: Brain mask dimensions do not match PET image dimensions.\n")
            return()
          }
          if (!all(dim(control_mask_nifti) == dim(first_pet_nifti))) {
            consoleMessages("Error: Control mask dimensions do not match PET image dimensions.\n")
            return()
          }
          
          # Create output directory
          base_output_dir <- file.path(getwd(), "Output_RAVEL_Normalization")
          dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)
          
          # Apply RAVEL normalization with error handling
          ravel_results <- tryCatch({
            normalizeRAVEL(
              input.files = pet_files,
              control.mask = control_mask_file,
              brain.mask = brain_mask_file,
              returnMatrix = TRUE,
              writeToDisk = FALSE,
              WhiteStripe = FALSE,
              verbose = TRUE
            )
          }, error = function(e) {
            consoleMessages(paste0("Error during RAVEL normalization: ", e$message, "\n"))
            return(NULL)
          })
          
          # Check if RAVEL returned valid results
          if (!is.null(ravel_results) && is.matrix(ravel_results)) {
            # Initialize lists to store input and output file paths
            input_files <- list()
            output_files <- list()
            
            # Process each PET image
            for (i in seq_along(pet_files)) {
              pet_file <- pet_files[i]
              
              # Read the input PET image as NIfTI
              pet_nifti <- readNIfTI(pet_file, reorient = FALSE)
              
              # Get the normalized intensity vector for this image
              normalized_intensity_vector <- ravel_results[, i]
              
              # Ensure the length of normalized_intensity_vector matches the number of voxels in brain mask
              num_mask_voxels <- sum(brain_mask_nifti)
              if (length(normalized_intensity_vector) != num_mask_voxels) {
                consoleMessages(paste0("Length of normalized intensities does not match number of masked voxels for ", basename(pet_file), ". Skipping.\n"))
                next
              }
              
              # Assign normalized intensities and set background to zero
              normalized_nifti <- pet_nifti
              normalized_nifti[brain_mask_nifti == 0] <- 0  # Set background to zero
              normalized_nifti[brain_mask_nifti == 1] <- normalized_intensity_vector
              
              # Save the normalized NIfTI image
              output_file <- file.path(base_output_dir, paste0("RAVEL_normalized_", sub("\\.nii\\.gz$", "", basename(pet_file))))
              
              
              writeNIfTI(normalized_nifti, filename = output_file)
              consoleMessages(paste0("Saved RAVEL normalized image: ", output_file, "\n"))
              
              # Store input and output file paths
              input_files[[basename(output_file)]] <- pet_file
              output_files[[basename(output_file)]] <- output_file
            }
            
            consoleMessages(paste0("RAVEL Normalization completed. All outputs saved to ", base_output_dir, "\n"))
            
            # --- RAVEL report (per output file) ---
            out_file_full <- paste0(output_file, ".nii.gz")  # writeNIfTI appends .nii.gz usually
            in_img  <- readNIfTI(pet_file, reorient = FALSE)
            out_img <- readNIfTI(out_file_full, reorient = FALSE)
            
            report_html <- generate_processing_report(
              step = "RAVEL Normalization",
              output_dir = base_output_dir,  # OR per-file folder (recommended)
              input_paths = list(
                PET = pet_file,
                BrainMask = brain_mask_file,
                ControlMask = control_mask_file
              ),
              output_paths = list(RAVEL_Normalized_PET = out_file_full),
              parameters = list(
                method = "normalizeRAVEL",
                WhiteStripe = FALSE,
                writeToDisk = FALSE,
                returnMatrix = TRUE
              ),
              logs = consoleMessages(),
              input_img = in_img,
              output_img = out_img
            )
            
            consoleMessages(paste0(consoleMessages(), "Report saved: ", report_html, "\n"))
            
            
            # Update the selectInput with processed outputs
            updateSelectInput(session, "selected_output",
                              choices = names(output_files),
                              selected = names(output_files)[1])
            
            # Initialize reactive values to store the images
            input_image <- reactiveVal()
            normalized_image <- reactiveVal()
            
            # Load the images based on the selected output
            observeEvent(input$selected_output, {
              selected_name <- input$selected_output
              req(selected_name)
              
              # Load the input and normalized images
              pet_file <- input_files[[selected_name]]
              output_file <- output_files[[selected_name]]
              
              input_image(readNIfTI(pet_file, reorient = FALSE))
              normalized_image(readNIfTI(output_file, reorient = FALSE))
              
              # Update slider input max value based on image dimensions
              img_dims <- dim(normalized_image())
              updateSliderInput(session, "slice_number", min = 1, max = img_dims[3],
                                value = round(img_dims[3] / 2))
            })
            
            # Load images for the first selection
            first_selection <- names(output_files)[1]
            if (!is.null(first_selection)) {
              selected_name <- first_selection
              pet_file <- input_files[[selected_name]]
              output_file <- output_files[[selected_name]]
              
              input_image(readNIfTI(pet_file, reorient = FALSE))
              normalized_image(readNIfTI(output_file, reorient = FALSE))
              
              # Update slider input max value based on image dimensions
              img_dims <- dim(normalized_image())
              updateSliderInput(session, "slice_number", min = 1, max = img_dims[3],
                                value = round(img_dims[3] / 2))
            }
            
            # Observers for navigating the slices
            observeEvent(input$prev_slice, {
              play_state(FALSE) # Stop play when manually changing slices
              updateActionButton(session, "play_slices", disabled = FALSE)
              updateActionButton(session, "stop_slices", disabled = TRUE)
              current_slice <- input$slice_number
              if (current_slice > 1) {
                updateSliderInput(session, "slice_number", value = current_slice - 1)
              }
            })
            
            observeEvent(input$next_slice, {
              play_state(FALSE) # Stop play when manually changing slices
              updateActionButton(session, "play_slices", disabled = FALSE)
              updateActionButton(session, "stop_slices", disabled = TRUE)
              current_slice <- input$slice_number
              if (current_slice < dim(normalized_image())[3]) {
                updateSliderInput(session, "slice_number", value = current_slice + 1)
              }
            })
            
            # Play slices
            observeEvent(input$play_slices, {
              play_state(TRUE)
              updateActionButton(session, "play_slices", disabled = TRUE)
              updateActionButton(session, "stop_slices", disabled = FALSE)
            })
            
            # Stop slices
            observeEvent(input$stop_slices, {
              play_state(FALSE)
              updateActionButton(session, "play_slices", disabled = FALSE)
              updateActionButton(session, "stop_slices", disabled = TRUE)
            })
            
            # Automatically update the slice number during play state
            observe({
              if (play_state()) {
                current_slice <- input$slice_number
                if (current_slice < dim(normalized_image())[3]) {
                  updateSliderInput(session, "slice_number", value = current_slice + 1)
                  invalidateLater(100, session)  # Adjust the interval as needed
                } else {
                  play_state(FALSE) # Stop play when the last slice is reached
                  updateActionButton(session, "play_slices", disabled = FALSE)
                  updateActionButton(session, "stop_slices", disabled = TRUE)
                }
              }
            })
            
            # Render the plot
            output$brain_plot <- renderPlot({
              req(input_image(), normalized_image())
              slice_num <- as.integer(input$slice_number)
              
              # Extract the slice as a matrix
              original_slice <- input_image()[,,slice_num]
              normalized_slice <- normalized_image()[,,slice_num]
              
              # Ensure slices are matrices
              if (!is.matrix(original_slice) || !is.matrix(normalized_slice)) {
                consoleMessages("Error: Extracted slices are not matrices.\n")
                return(NULL)
              }
              
              # Set up the plotting area to have 1 row and 2 columns
              par(mfrow = c(1, 2), oma = c(2, 1, 2, 2), mar = c(2, 2, 2, 2),
                  bg = "black", fg = "white")
              
              # Function to plot a single slice with black background if empty
              plot_slice <- function(slice, title) {
                if (all(slice == 0)) {
                  par(bg = "black")
                  plot.new()
                  plot.window(xlim = c(0, 1), ylim = c(0, 1))
                  rect(0, 0, 1, 1, col = "black", border = "black")
                  title(title, col.main = "white")
                } else {
                  image(t(apply(slice, 2, rev)), col = gray(0:64/64), axes = FALSE, main = "", asp = 1)
                  title(title, col.main = "white")
                }
              }
              
              # Display the selected slice for the input and normalized images
              plot_slice(original_slice, "Original PET Image")
              plot_slice(normalized_slice, "RAVEL Normalized Image")
              
              # Add filename at the bottom
              #selected_name <- input$selected_output
              # mtext(paste("Input Filename:", selected_name), side = 1, line = -2, adj = 0.5, cex = 1.5, col = "white", outer = TRUE)
            }, res = 200)
            # Reset all inputs after brain extraction is complete
            updateTextInput(session, "pet_input", value = "")
            updateTextInput(session, "mri_input", value = "")
            updateTextInput(session, "csv_input", value = "")
            updateTextInput(session, "control_mask_input", value = "")
            updateTextInput(session, "brain_mask_input", value = "")
            
          } else {
            consoleMessages("Error: RAVEL normalization failed. Please check the input files and masks.\n")
          }
        } else {
          consoleMessages(paste0(consoleMessages(), "Please provide PET scans, control mask, and brain mask for RAVEL Normalization.\n"))
        }
      }
      
      
      incProgress(0.98, detail = "Done")
    })
    
  })
  
  output$console <- renderText({
    consoleMessages()
  })
  
  # Update text input with chosen file for PET Scans Input Path
  observeEvent(input$btn_pet, {
    if (!is.null(input$btn_pet)) {
      pet_path <- shinyFiles::parseFilePaths(volumes, input$btn_pet)
      if (nrow(pet_path) > 0) {
        # Join multiple file paths or just show the folder path
        pet_paths <- paste(pet_path$datapath, collapse = ",")
        updateTextInput(session, "pet_input", value = pet_paths)
        consoleMessages(paste0(consoleMessages(), "PET scan(s) selected: ", dirname(pet_path$datapath[1]), "\n"))
      }
    }
  })
  
  # Update text input with chosen file for MRI Scans Input Path
  observeEvent(input$btn_mri, {
    if (!is.null(input$btn_mri)) {
      mri_path <- shinyFiles::parseFilePaths(volumes, input$btn_mri)
      if (nrow(mri_path) > 0) {
        # Join multiple file paths or just show the folder path
        mri_paths <- paste(mri_path$datapath, collapse = ",")
        updateTextInput(session, "mri_input", value = mri_paths)
        consoleMessages(paste0(consoleMessages(), "MRI scan(s) selected: ", dirname(mri_path$datapath[1]), "\n"))
      }
    }
  })
  
  # Update text input with chosen file for CSV Input Path
  observe({
    if (!is.null(input$btn_csv)) {
      csv_path <- shinyFiles::parseFilePaths(volumes, input$btn_csv)
      if (nrow(csv_path) > 0) {
        # Convert the path to a character string
        csv_file_path <- as.character(csv_path$datapath[1])
        # Update the text input with the string path
        updateTextInput(session, "csv_input", value = csv_file_path)
      }
    }
  })
}

# Function to get an available port
get_available_port <- function() {
  port <- httpuv::randomPort()
  return(port)
}

# Get an available port
port <- get_available_port()
message("Selected port:", port)

# Run the Shiny app
shinyApp(
  ui = ui,
  server = server
)
