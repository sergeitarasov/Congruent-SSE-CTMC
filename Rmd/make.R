
bookdown::render_book("Rmd", output_format = "bookdown::gitbook", output_dir = '../Report')

#bookdown::render_book("Rmd", output_format = "bookdown::pdf_book", output_dir = '../Report-pdf')

I=matrix(c(1,2,3,4), 2,2, byrow = T)
Fi= I=matrix(c(2,0,0,2), 2,2, byrow = T)
I %*% Fi
