(TeX-add-style-hook
 "rr"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "a4paper")))
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

