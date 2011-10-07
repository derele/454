## Common R function fro 454-analysis and beyond
## Emanuel Heiltlinger: emanuelheitlinger@gmail.com

read.sequences <- function (file){
  content <- readLines(file)
  header.lines <- grep( "^>", content)
  start.lines <- header.lines+1
  end.lines <- c(header.lines[-1]-1, length(content))
  sq <- sapply(1:length(start.lines), function (i) {
    list(content[start.lines[i]:end.lines[i]])})
  names(sq) <- substr(content[header.lines],2, nchar(content[header.lines]))
  sq <- unlist(lapply(sq, paste, collapse=""))
  return(sq[nchar(sq)>0])
}

write.sequence <- function (sequence.obj, path){
  write(paste(">",
              names(sequence.obj), "\n",
              sequence.obj, sep=""),
        path)
}

read.blast.best <- function (path){ 
  command <- paste ("rev", path, "| uniq -f 11 | rev")
  read.delim(pipe(command), header=FALSE)}

reduce.blast.by.length <- function (bt, len=0.8, iden=95){
br <- bt[(bt$length*len < bt$e.hit-bt$s.hit & bt$iden>iden), ]
br <- br[order(br[,1], br[,5], decreasing=TRUE),]
br <- br[!duplicated(br[,1]),]
}


get.gc <- function (x) (nchar(gsub("A|T", "",  x))/nchar(as.character(x)))*100

strReverse <- function(y) sapply(lapply(strsplit(y, NULL), rev), paste,
                                 collapse="")

revcom <- function(x){
  r <- strReverse(x)
  return(chartr("acgtryswkmbdhvnxACGTRYSWKMBDHVNX", "tgcayrswmkvhdbnxTGCAYRSWMKVHDBNX", r))
}



