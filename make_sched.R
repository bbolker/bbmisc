library(googledrive)
library(gmailr)
## unnecessary (maybe useful in the future for fine control?
## library(googlesheets4)

## TODO
## * clean up / modularize / document
## * store output as tempfile (for safety)?
## * 'update' option?

prefix <- "sched_winter2020"
people_fn <- sprintf("%s_people.csv",prefix)   ## people to include (name, email)
exclude_fn <- sprintf("%s_exclude.csv",prefix) ## slots to exclude (day, time, reason)
out_fn <- paste0(prefix,".csv")

msg <- "Please enter 'yes' (default), 'no', or 'maybe' to indicate your availability during a typical week Jan-May 2020. (All times are America/Toronto.)"
## cc <- read.csv("classlist.csv", stringsAsFactors=FALSE)
## email <- gsub("@.*$","",cc$Email)

options(stringsAsFactors=FALSE)
weekdays <-  c("Mon","Tues","Weds","Thurs","Fri")

##' @param hour_start minutes for hour start
##' @param days character
##' @param hours numeric
expand_slots <- function(hour_start=30,
                         days=weekdays,
                         hours=c(9:12,1:4)) {
    slots <- sprintf("%s:%d",hours,hour_start)
    return(apply(expand.grid(slots=slots,days=days)[c("days","slots")],
                 1,paste,collapse=" "))
}


#' @param time slots
#' @param IDs people
#' @param default_val default cell entry
#' @param vert_par which parameter goes top-to-bottom?
make_sched0 <- function(slots, IDs, default_val="",
                        vert_par=c("times","IDs")) {
    vert_par <- match.arg(vert_par)
    mm <- matrix(default_val,
                 nrow=length(slots),ncol=length(IDs),
                 dimnames=list(slot=slots,ID=IDs))
    if (vert_par=="IDs") mm <- t(mm)
    return(as.data.frame(mm))
}

make_sched <- function(IDs="") {
    make_sched0(slots=expand_slots(),IDs=IDs)
}

#' post-process to add information about excluded times
exclude_times <- function(s,exclude_data, action=c("label","delete")) {
    action <- match.arg(action)
    if (action=="delete") stop("'delete' option not yet implemented")
    times <- with(exclude_data,paste(day,time))
    for (i in seq(nrow(exclude_data))) {
        s[times[i],] <- sprintf("no\n(%s)",exclude_data$reason[i])
    }
    return(s)
}
    
## (obsolete) stuff for Britton lecture
## ee <- expand_slots(days=paste(c("Tues","Weds","Thurs","Fri"),
##                               "Apr",9:12),
##                    hours=c(9:12,1,6))
## ee <- gsub("6:30","dinner",ee)
## ee <- gsub("12:30","12:30 (lunch)",ee)
## ss_britton <- make_sched0(slots=ee,IDs="Ben Bolker")
## write.csv(ss_britton,file="sched_britton.csv")
## sched_ss <- gs_upload("sched_britton.csv")

people <- read.csv(people_fn)
exclude <- read.csv(exclude_fn)

ss <- make_sched(people$name)
ss <- exclude_times(ss,exclude)
write.csv(ss,file=out_fn)
View(ss)

sched_ss <- drive_upload(out_fn)


## uploaded into Drive file:
##   * sched_winter2020.csv: 14XDgjWRthwB1usWbRwi-HqAVoQzDgDnO

## https://cran.r-project.org/web/packages/googlesheets/vignettes/basic-usage.html#make-new-sheets-from-local-delimited-files-or-excel-workbooks

lnk <- drive_link(sched_ss)

for (email in people$email) {
    sched_ss <- drive_share(sched_ss,
                            role = "writer",
                            type = "user",
                            emailAddress = email,
                            emailMessage = msg
         )
}

browseURL(lnk)
