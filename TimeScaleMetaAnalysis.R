library(shiny)
library(expm)
library(nloptr)
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),  # Load the shinyjs library
  titlePanel("Matrix and Delta Input"),
  sidebarLayout(
    sidebarPanel(
      numericInput("num_matrices", "Number of Matrices:", value = 1, min = 1),
      numericInput("num_dimensions", "Number of Dimensions:", value = 2, min = 1, max = 10),
      actionButton("confirm_button", "Confirm"),
      uiOutput("matrix_input")
    ),
    mainPanel(
      textOutput("confirmation_text"),
      div(style="display: flex; justify-content: center;",
          tableOutput("final_matrix_table")  # Center the table
      ),
      plotOutput("matrix_plot")
    )
  )
)

server <- function(input, output, session) {
  num_matrices <- reactiveVal()
  num_dimensions <- reactiveVal()
  current_matrix <- reactiveVal(1)
  
  observeEvent(input$confirm_button, {
    num_matrices(input$num_matrices)
    num_dimensions(input$num_dimensions)
    current_matrix(1)
    
    shinyjs::disable("num_dimensions")  # Disable the numeric input
  })
  
  output$matrix_input <- renderUI({
    if (!is.null(num_matrices()) && current_matrix() <= num_matrices()) {
      matrix_dim <- numericInput(inputId = paste0("matrix_dim_", current_matrix()), 
                                 label = paste("Matrix", current_matrix(), "Dimension (Rows):"), 
                                 value = num_dimensions(), min = 1, max = 10)
      
      matrix_values <- textInput(inputId = paste0("matrix_values_", current_matrix()), 
                                 label = paste("Matrix", current_matrix(), "Values (Comma-Separated):"), 
                                 value = "")
      
      matrix_delta <- numericInput(inputId = paste0("matrix_delta_", current_matrix()), 
                                   label = paste("Delta for Matrix", current_matrix()), value = 1)
      
      action_button <- actionButton("next_button", "Next Matrix")
      
      matrix_inputs <- tagList(matrix_dim, matrix_values, matrix_delta, action_button)
      return(matrix_inputs)
    }
  })
  
  observeEvent(input$next_button, {
    if (current_matrix() < num_matrices()) {
      current_matrix(current_matrix() + 1)
    } else {
      output$confirmation_text <- renderText({
        "Matrices and delta values confirmed."
      })
      
      # Create a list of matrices and deltas based on user input
      matrices <- list()
      deltas <- list()
      NE <- 0
      
      for (i in seq_len(num_matrices())) {
        matrix_dim <- num_dimensions()
        matrix_values <- as.numeric(strsplit(input[[paste0("matrix_values_", i)]], ",")[[1]])
        matrix_values <- matrix(matrix_values, nrow = matrix_dim, byrow = TRUE)
        matrices[[paste0("mat", i)]] <- matrix_values
        deltas[[paste0("del", i)]] <- input[[paste0("matrix_delta_", i)]]
        NE <- matrix_dim^2
      }
      
      # Initial values for optimization
      params <- rnorm(NE, 0, 0.1)
      
      # Objective Function to be Minimized
      eval_f0 <- function(params, matrices, deltas) {
        # Reshape params into a matrix
        target_matrix <- matrix(params, nrow = sqrt(NE), ncol = sqrt(NE), byrow = TRUE)
        # Compute the sum of absolute differences between matrices
        distances = matrix(0, length(matrices), 1)
        for(i in 1:length(matrices)){
          distances[i,] = sum(abs(expm(target_matrix*deltas[[i]]) - matrices[[i]]))
        }
        distance = sum(distances)
        return(distance)
      }
      
      # Options for nloptr
      opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD",
                   "xtol_rel" = 1e-6,
                   "maxeval" = 2000)
      
      # Running the Optimization
      output_optimization <- nloptr(x0 = params, eval_f = eval_f0,
                                     matrices = matrices, deltas = deltas,
                                     opts = opts,
                                     lb = rep(-0.9, NE),
                                     ub = rep(0.9, NE))
      
      # Construct the final matrix
      final_matrix <- matrix(output_optimization$solution, nrow = sqrt(NE), ncol=sqrt(NE), byrow = TRUE)
      
      # Set row and column names
      rownames(final_matrix) = colnames(final_matrix) = paste0('V', 1:sqrt(NE))
      
      # Display the final matrix as an HTML table
      output$final_matrix_table <- renderTable({
        final_matrix
      })
      
      # Update the plot with the optimized matrix and legend
      output$matrix_plot <- renderPlot({
        # Create plot_mat
        plot_mat <- matrix(final_matrix, nrow = sqrt(NE), ncol = sqrt(NE), byrow = TRUE)
        
        # Plotting OU by Delta-T
        delts = seq(0.0, 20, 0.5)
        transforms = matrix(NA, length(delts), NE)
        index = 1
        for(i in delts){
          transforms[index,] = as.numeric(expm(plot_mat*i))
          index = index + 1
        }
        
        paramnames = c(matrix(levels(interaction(letters[1:sqrt(NE)], 1:sqrt(NE), sep = "_")),
                              sqrt(NE), sqrt(NE),byrow = TRUE))
        
        plot(delts, transforms[,1], type = 'n', ylim = c(0,1.1), lwd = 1.5, lty = 2, 
             ylab = 'Coefficient Value', xlab = expression(Delta * t))
        for (i in 1:NE) {
          lines(delts, transforms[,i], type = 'l', col = rainbow(NE)[i], lwd = 1.5)
        }
        
        legend("topright", legend = paste(paramnames), col = rainbow(NE), lwd = 1.5, cex = 0.8)
      })
    }
  })
}

shinyApp(ui, server)
