library(ggplot2); theme_set(theme_bw())
today <- Sys.Date()

sscrape <- function(string="logistic+regression") {
    require("stringr")
    sstring0 <- "http://scholar.google.ca/scholar?as_q=&num=10&btnG=Search+Scholar&as_epq=STRING&as_oq=&as_eq=&as_occt=any&as_sauthors=&as_publication=&as_ylo=&as_yhi=&as_sdt=1.&as_sdtp=on&as_sdts=5&hl=en"
    sstring <- sub("STRING",string,sstring0)
    rr <- suppressWarnings(readLines(url(sstring), encoding = "utf8"))
    ## rr2 <- rr[grep("[Rr]esults",rr)[1]]
    ## rr2 <- head(rr, 10) ## don't need the whole thing, avoid encoding problems
    ## rstr <- rr2 |>
    ##     gsub(pattern = "^.+[Rr]esults.+of about <b>", replacement = "") |>
    ##     gsub(pattern = ",", replacement = "") |>
    ##     gsub(pattern = "</b>.+$", replacement = "")
    rstr <- na.omit(stringr::str_extract(rr,"About [0-9,]+ results"))
    rnum <- as.numeric(gsub(",","",str_extract(rstr,"[0-9,]+")))
    attr(rnum,"scrape_time") <- Sys.time()
    return(rnum)
}


fn <- "gscrape.RData"
## could use a caching solution for Sweave (cacheSweave, weaver package,
##  pgfSweave ... but they're all slightly wonky with keep.source at
##  the moment
search_terms <- c("generalized+linear+model",
                  "logistic+regression",
                  "Poisson+regression",
                  "binomial+regression")
if (!file.exists(fn)) {
  gscrape <- sapply(search_terms, sscrape)
  save("gscrape",file=fn)
} else load(fn)

d <- data.frame(n=names(gscrape),v=gscrape)
d$n <- reorder(d$n,d$v)
d <- d[order(d$v),]
gg1 <- ggplot(d,aes(x=v,y=n))+geom_point(size=5)+
    ## xlim(0.5e4,2e6)+
    scale_x_log10() +
    geom_text(aes(label=v),colour="red",vjust=2)+
    labs(y="",x="Google Scholar hits", subtitle = paste("Search conducted on", today))

alt_text <- with(d, sprintf("%s: %d hits", n, v)) |> paste(collapse = "; ")
alt_text <- paste0("Results of a Google Scholar search on ",
                   format(today),": ", alt_text)

print(gg1)

ggsave("glm_scrape.png")
