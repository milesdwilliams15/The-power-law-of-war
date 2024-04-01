# combine the title page and manuscrip
# pdftools::pdf_combine(
#   c(here::here("04_report", "01_title_page.pdf"),
#     here::here("04_report", "01_manuscript.pdf")),
#   here::here("04_report", "01_combined_manuscript.pdf")
# )
pdftools::pdf_subset(
  input = here::here("04_report", "04_manuscript.pdf"),
  output = here::here("04_report", "04_title.pdf"),
  pages = 1
)
pdftools::pdf_subset(
  input = here::here("04_report", "04_manuscript.pdf"),
  output = here::here("04_report", "04_anon_man.pdf"),
  pages = 2:54
)
