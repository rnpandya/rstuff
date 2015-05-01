library(synapseClient)
synapseLogin()

synBigQuery <- function (query) {
	limit <- 1000
	while (limit > 0) {
		offset = 1
		tryCatch({
			result <- NULL
			while (T) {
				cat(paste(query, "limit", limit, "offset", offset, "\n"))
				chunk <- synQuery(paste(query, "limit", limit, "offset", offset))
				result <- if (is.null(result)) chunk else rbind(result, chunk)
				if (dim(chunk)[1] < limit) {
					return(result)
				}
				offset <- offset + limit
			}
		}, error = function (e) {
			cat(paste("error:", as.character(e), "\n"))
		})
		limit <- as.integer(limit / 10)
	}
	stop("unable to run query", query)
}
