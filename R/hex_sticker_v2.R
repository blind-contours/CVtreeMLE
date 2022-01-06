#' Creates the CVtreeMLE sticker from the image in figures
#'
#' @import hexSticker
#' @import here
#'
#' @export
col_bg <- "#f5ab35"      ## Lighning yellow
col_border <- "#f5d76e"  ## Cream can
col_text <- "#22313f"    ## Ebony clay

hexSticker::sticker(

  # image
  here::here("inst/figures/CVtreeMLE.001.png"),
  s_x=1.00, # slightly to right to appear centered
  s_y=0.98,
  s_width=0.98,
  s_height=.99,

  # package name
  package="CVtreeMLE",
  p_family = "Aller_Lt",
  p_size=6.1,
  p_color = col_text,
  h_fill = col_bg,
  h_color = col_border,
  p_x = 1.0,
  p_y = 1.51,

  # Output file
  filename=here::here("inst/figures/CVtreeMLE_sticker.png"),

  # Background colour

  # Border
  # Grey colours: https://www.w3schools.com/colors/colors_shades.asp
  # h_size = 1.5,

  dpi = 2000 # otherwise the final fantasy image quality is not good
)
