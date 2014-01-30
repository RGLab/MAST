\name{readNanoStringLanes}
\alias{readNanoStringLanes}
\title{RCC Nanostring Reader}
\usage{
readNanoStringLanes(x)
}
\arguments{
  \item{x}{is a \code{character} vector of full path names
  to RCC files}
}
\value{
a list
}
\description{
Reads RCC Nanostring Lane files and constructs a
SingleCellAssay object
}
\details{
Open each lane file and read line by line constructing the
different components of the file

Metadata is present in the files. One lane per file.
}

