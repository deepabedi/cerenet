# Suppress "no visible binding" notes for these variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("imagecol", "imagerow", ".data"))
}

