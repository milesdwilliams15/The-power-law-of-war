# combine the title page and manuscrip
pdftools::pdf_combine(
  c(here::here("04_report", "01_title_page.pdf"),
    here::here("04_report", "01_manuscript.pdf")),
  here::here("04_report", "01_combined_manuscript.pdf")
)
