(TeX-add-style-hook
 "zeta_ophiuchi"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("aastex63" "twocolumn" "twocolappendix" "trackchanges")))
   (TeX-run-style-hooks
    "latex2e"
    "aastex63"
    "aastex6310"
    "amsmath"
    "CJK")
   (TeX-add-symbols
    '("todo" 1)
    '("Secref" 1)
    '("Tabref" 1)
    '("Figref" 1)
    '("Eqref" 1)
    '("code" 1)
    "mesa"
    "MESA"
    "kms"
    "Msun"
    "kev"
    "gk")
   (LaTeX-add-labels
    "sec:intro"
    "sec:methods"
    "sec:software")
   (LaTeX-add-bibliographies
    "./zeta_oph.bib"))
 :latex)

