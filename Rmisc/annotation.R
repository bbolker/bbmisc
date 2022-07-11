## adapts ggimage::geom_subview to pass data (facet info)
geom_subview <- function (mapping = NULL, data = NULL, width = 0.1, height = 0.1, 
    x = NULL, y = NULL, subview = NULL, panel = NULL)
{
    get_aes_var <- ggfun::get_aes_var
    if (is.null(data)) {
        data <- tibble(x = x, y = y)
    }
    else if (!inherits(data, "tbl")) {
        data <- as_tibble(data)
    }
    if (is.null(mapping)) {
        mapping <- aes_(x = ~x, y = ~y)
    }
    mapping <- as.list(mapping)
    if (is.null(mapping$x)) {
        stop("x aesthetic mapping should be provided")
    }
    if (is.null(mapping$y)) {
        stop("y aesthetic mapping should be provided")
    }
    if (is.null(mapping$subview) && is.null(subview)) {
        stop("subview must be provided")
    }
    if (is.null(mapping$subview)) {
        if (!inherits(subview, "list")) {
            subview <- list(subview)
        }
        data$subview <- subview
    }
    else {
        sv_var <- get_aes_var(mapping, "subview")
        data$subview <- data[[sv_var]]
    }
    xvar <- get_aes_var(mapping, "x")
    yvar <- get_aes_var(mapping, "y")
    if (is.null(mapping$width)) {
        data$width <- width
    }
    else {
        width_var <- get_aes_var(mapping, "width")
        data$width <- data[[width_var]]
    }
    if (is.null(mapping$height)) {
        data$height <- height
    }
    else {
        height_var <- get_aes_var(mapping, "height")
        data$height <- data[[height_var]]
    }
    data$xmin <- data[[xvar]] - data$width/2
    data$xmax <- data[[xvar]] + data$width/2
    data$ymin <- data[[yvar]] - data$height/2
    data$ymax <- data[[yvar]] + data$height/2
    lapply(1:nrow(data), function(i) {
        annotation_custom(as.grob(data$subview[[i]]), xmin = data$xmin[i], 
                          xmax = data$xmax[i], ymin = data$ymin[i], ymax = data$ymax[i],
                          data = data[i,])
    })
}

## annotation_custom
## includes a data argument so we can pass facet/panel data in
annotation_custom <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                               ##panel = NULL
                               data = ggplot2:::dummy_data()
                               )
{
    layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = GeomCustomAnn, inherit.aes = FALSE, params = list(grob = grob, 
            xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax))
}

