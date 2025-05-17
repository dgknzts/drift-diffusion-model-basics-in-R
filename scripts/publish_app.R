# Script to publish Shiny app to shinyapps.io

# Install required packages if needed
if (!require("rsconnect")) {
  install.packages("rsconnect")
  library(rsconnect)
}

# Configuration - only needed first time
# Uncomment and fill in these lines the first time you run this script
rsconnect::setAccountInfo(name='dgknzts',
			  token='97D31B701FF0FFCF6F69F0AF1016C417',
			  secret='PpiCMJoComA6xT9+P7dMAqu2rs7QCqz9qzc9RBz5')

# Publish the app
rsconnect::deployApp(
  appDir = ".",           # Current directory with app.R
  appName = "DDM-Basics",  # Name for your app on shinyapps.io
  appTitle = "DDM Basics Interactive Demo"
)

cat("App published to shinyapps.io!\n")
cat("After publishing, you'll get a URL that looks like:\n")
cat("https://dgknzts.shinyapps.io/DDM-Basics/\n")
cat("Use this URL to add to your pkgdown site.\n") 