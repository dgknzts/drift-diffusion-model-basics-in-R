# Script to update pkgdown site when new content is added

# Clear previous build artifacts (optional but recommended)
if (dir.exists("docs/articles")) {
  unlink("docs/articles/*.html")
}

# Build the pkgdown site
pkgdown::build_articles()

# Print message when done
cat("PKGdown site updated successfully!\n")
cat("Remember to commit and push changes to GitHub:\n")
cat("git add .\n")
cat("git commit -m \"Update site with new vignette\"\n")
cat("git push\n") 