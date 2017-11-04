(TeX-add-style-hook
 "rr"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "packages"
    "macro"
    "body"
    "article"
    "art10"
    "RR"
    "a4wide"
    "natbib")
   (LaTeX-add-bibliographies
    "./biblio/String"
    "./biblio/NonSmooth"
    "./biblio/Math"
    "./biblio/Multibody"
    "./biblio/Fem.bib"
    "./biblio/Dae.bib"
    "./biblio/Meca"
    "./biblio/AnaNum.bib"
    "./biblio/Math-Impact"
    "./biblio/Contact"
    "./biblio/Optim"
    "./biblio/Cp"))
 :latex)

