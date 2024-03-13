## https://stackoverflow.com/questions/50016207/unable-to-send-email-using-mailr-package

## https://myaccount.google.com/lesssecureapps

## run from head of private repo
## need to enable "insecure mode" in Gmail
## Google Account -> Security

## app passwords

## https://support.google.com/accounts/answer/185833?hl=en
library(dplyr)
## library(googlesheets4)
library(mailR)
library(getPass)
library(readr)

## Google sheets code
## gg <- gs_title("Math3MB3_marks ")
## gg <- gs_key("1q9jivlvZ0MREyGdw8quyxibjCJjLmjrABsAyye7xjrw")
## gg2 <- gs_read(gg,ws="project",check.names=TRUE) %>%
##     select(Last.Name,First.Name,Email,topic,group_2) %>%
##     filter(!is.na(group_2))
## submdir <- "."
## subm <- list.files(submdir,pattern="\\.txt$")


## grps <- gg2 %>% group_by(group_2) %>%
##     mutate(fullname=paste(First.Name,Last.Name)) %>%
##     summarise(all_partic=paste(fullname,collapse=", "))

## body <- readLines("init_proposal.txt")

send_mail <- function(subject,
                      recipient,
                      body_text,
                      body_file,
                      fake=TRUE,
                      attach_files,
                      gpass = NULL,
                      sender="bbolker@gmail.com",
                      user.name=gsub("@.*$","",sender),
                      port=465,
                      ssl=TRUE,
                      host.name="smtp.gmail.com"
                      ) {

    if (!missing(body_text)) {
        body_file <- "body.tmp"
        writeLines(body_text,con=body_file)
    }

  ## cache password **in global working dir** ... sloppy ...
  if (is.null(gpass)) {
    if (!is.null(.gpass))
      gpass <- .gpass
  }
  if (gpass == "RESET" || is.null(gpass)) {
    gpass <- getPass::getPass()
    .gpass <<- gpass
  }

  ## FIXME: multiple attachments?
  if (!is.null(attach_files) && !is.na(attach_files)) {
    attach_files <- trimws(strsplit(attach_files, ";")[[1]])
    missing <- which(!sapply(attach_files, file.exists))
    if (length(missing)>0) {
      stop("missing attachments: ", paste(names(missing), collapse=", "))
    }
  }

  if (length(attach_files) == 1 && is.na(attach_files)) {
    attach_files <- NULL
  }

  mailR::send.mail(from = sender,
              to = recipient,
              subject = subject,
              body = body_file,
              smtp = list(host.name = host.name,
                          port = port,
                          user.name = user.name,
                          passwd = gpass,
                          ssl = ssl),
              authenticate = TRUE,
              attach.files = attach_files,
              send = !fake)
    cat(recipient,"\n")
}

.gpass <- NULL

## sub-tasks:
##   find files with appropriate tags and extensions in subdir
##   substitute strings in body text

substitute_strings <- function(text, replacements) {
  replacements <- as.list(replacements)
  for (i in seq_along(replacements)) {
    text <- gsub(sprintf("\\$\\{%s\\}", names(replacements)[i]), replacements[[i]], text)
  }
  return(text)
}

## vectorized search for a single file matching each pattern in a string
##' @param pat pattern
##' @param exclude pattern to exclude
##' @param allow_mult allow multiple files per individual?
##' @param dir directory to search
get_files <- function(pat, exclude = NULL,
                      allow_mult = c("error", "latest", "OK"), dir = ".") {
    allow_mult <- match.arg(allow_mult)
    res <- vector("list", length(pat))
    names(res) <- pat
    for (p in pat) {
        f <- list.files(path=dir, pattern = p, ignore.case = TRUE)
        if (length(exclude) > 0) {
            f <- f[!grepl(paste(exclude, collapse="|"), f)]
        }
        if (length(f) > 1) {
            if (allow_mult == "error") stop("multiple files: ",paste(f,collapse=", "))
            if (allow_mult == "latest") {
                f <- f[which.max(sapply(f, function(x) file.info(x)[["mtime"]]))]
            }
        }
        f <- paste(f, collapse = "; ")
        res[[p]] <- if (length(f)==0) NA else f
    }
    return(res)
}



#' @param docdir directory to search for documents
#' @param data info on mailing lists etc.
#' @param ... info to pass to send_mail (subject, ...)
#' @param do_self_test
#' @param body_file
#' @param body_text
do_mail <- function(data,
                    doc_dir = NULL,
                    doc_name = NULL,
                    self_test = TRUE,
                    self_email = "bbolker@gmail.com",
                    email_name = "Email",
                    mail_args = NULL,
                    skip_nodoc = FALSE,
                    body_text = NULL,
                    fake_run = TRUE,
                    ...) {
  if (!exists("gpass")) gpass <- getPass::getPass()

  get_attach <- function(doc_name, data) {
    if (is.null(doc_name) || is.na(data[[doc_name]]) || is.null(data[[doc_name]])) return(NULL)
    if (is.null(doc_dir)) return(data[[doc_name]])
    ## ugh!
    tmp <- sapply(strsplit(data[[doc_name]],";")[[1]],
                 function(x) file.path(doc_dir, x))
    paste(tmp, collapse = ";")
  }

  ## testing
  if (self_test) {
    first_data <- data[1,]
    first_doc_data <- if (is.null(doc_name)) NULL else data[which(!is.na(data[[doc_name]]))[1],]
    do.call(send_mail,
            c(list(...),
              list(
                  sender = self_email,
                  recipient = self_email,
                  body_text = substitute_strings(body_text, first_data),
                  gpass = gpass,
                  attach_file = get_attach(doc_name, first_doc_data),
                  fake = FALSE),
              mail_args))
  }

if (is.null(data[[email_name]])) stop("can't find email variable ", email_name)
  for (i in seq(nrow(data))) {
      if (!is.null(doc_name) && !skip_nodoc && is.na(fn <- data[[doc_name]][i])) next
      email <- data[[email_name]][i]
      if (is.na(email)) {
          cat("skipping entry ",i,": NA e-mail address\n")
          next
      }
      do.call(send_mail,
            c(list(...),
              list(
                   sender = self_email,
                   recipient = data[[email_name]][i],
                   body_text = substitute_strings(body_text, data[i,]),
                   gpass = gpass,
                   attach_file = get_attach(doc_name, data[i,]),
                   fake = fake_run),
              mail_args)
            )
  } ## loop over rows of data
}
