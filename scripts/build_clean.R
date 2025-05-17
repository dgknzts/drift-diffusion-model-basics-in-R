# Clean and rebuild pkgdown site
library(pkgdown)

# Clean existing files
if (dir.exists("docs/articles")) {
  file_list <- list.files("docs/articles", pattern = "\\.html$", full.names = TRUE)
  unlink(file_list)
  message("Cleaned existing HTML files in docs/articles/")
}

# Make sure any problem files in vignettes are removed
problem_file <- "vignettes/--find-assets.html"
if (file.exists(problem_file)) {
  unlink(problem_file)
  message("Removed problematic file: ", problem_file)
}

# Build the site with more control
message("Building articles...")
tryCatch({
  pkgdown::build_articles(quiet = FALSE)
  message("Articles built successfully!")
}, error = function(e) {
  message("Error building articles: ", e$message)
})

message("Done!") 