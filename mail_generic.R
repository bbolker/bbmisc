## https://stackoverflow.com/questions/50016207/unable-to-send-email-using-mailr-package

## https://myaccount.google.com/lesssecureapps

## run from head of private repo
## need to enable "insecure mode" in Gmail
## Google Account -> Security
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
                      attach_file,
                      gpass=getPass::getPass(),
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

    ## FIXME: multiple attachments?
    astr <- if (!is.null(attach_file) && file.exists(attach_file)) {
                attach_file
            } else NULL
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
              attach.files=astr,
              send = !fake)
    cat(recipient,"\n")
}

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
get_files <- Vectorize(
    vectorize.args = "pat",
    FUN =
      function(pat, exclude = NULL, allow_mult = FALSE) {
      f <- list.files(path=hwdir, pattern = pat, ignore.case = TRUE)
      if (length(exclude) > 0) {
        f <- f[!grepl(paste(exclude, collapse="|"), f)]
      }
      if (!allow_mult && length(f)>1) stop("multiple files: ",paste(f,collapse=", "))
      f <- paste(f, collapse = "; ")
      if (length(f)==0) return(NA)
      names(f) <- pat
      return(f)
    })


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
    file.path(doc_dir, data[[doc_name]])
  }

  ## testing
  if (self_test) {
    first_data <- data[1,]
    first_doc <- if (is.null(doc_name)) NULL else na.omit(data[[doc_name]])[1]
    do.call(send_mail,
            c(list(...),
              list(
                  sender = self_email,
                  recipient = self_email,
                  body_text = substitute_strings(body_text, first_data),
                  gpass = gpass,
                  attach_file = get_attach(first_doc, data),
                  fake = FALSE),
              mail_args))
  }

  for (i in seq(nrow(data))) {
    if (!is.null(doc_name) && !skip_nodoc && !is.na(fn <- data[[doc_name]][i])) next
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
