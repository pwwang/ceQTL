options(stringsAsFactors = FALSE)
library(chunked)
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
method  = {{args.fdr | quote}}
pval    = {{args.pval | R}}

# Case    Pooled  Groups  Fstat   Pval
# rs10774580.ZNF219.SPPL3 ZNF219=-0.257,_=8.809,N=703 AA(hom):ZNF219=-0.163,_=8.422,N=115; GA(het):ZNF219=-0.404,_=9.406,N=332; GG(ref):ZNF219=-0.078,_=8.089,N=256   4.754   8.650E-04
# rs3782410.TFAP2C.ITGA5  TFAP2C=-0.25,_=7.707,N=702  AA(hom):TFAP2C=-0.387,_=7.972,N=72; GA(het):TFAP2C=-0.09,_=7.305,N=302; GG(ref):TFAP2C=-0.372,_=8.032,N=328 5.047   5.159E-04
pvalues  = read.table(infile, header = TRUE, row.names = NULL, sep = "\t", colClasses = c("NULL", "NULL", "NULL", "NULL", "double"))
adjpvals = p.adjust(unlist(pvalues), method = method)

logger = function(...) {
	msg = paste(list(...), collapse = ' ')
	cat(paste0(msg, "\n"), file = stderr())
}

# in case in file is too large, read it in chunks
nlines    = length(adjpvals)
chunksize = 10000
logger('Chunk size:', chunksize)
nchunks   = ceiling(nlines / chunksize)
logger('No. of chunks:', nchunks)
chunks = read_chunkwise(infile,  chunk_size = chunksize, format = "table", header = TRUE)
sapply(1:nchunks, function(n) {
	logger('- Handling chunk #', n, '/', nchunks, '(', sprintf('%.2f%%', 100*n/nchunks), ') ...')
	start     = (n-1) * chunksize + 1
	chunkdata = chunks$next_chunk()
	outdata   = cbind(chunkdata, AdjPval = sprintf("%.3E", adjpvals[start:start + chunksize - 1]))
	outdata   = outdata[outdata$Pval < pval,,drop = FALSE]
	cnames    = n == 1
	write.table(outdata, outfile, append = !cnames, col.names = cnames, row.names = FALSE, sep = "\t", quote = FALSE)
})

