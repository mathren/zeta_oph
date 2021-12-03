(TeX-add-style-hook
 "zeta_ophiuchi"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("aastex63" "twocolumn" "twocolappendix" "trackchanges")))
   (TeX-run-style-hooks
    "latex2e"
    "table"
    "aastex63"
    "aastex6310"
    "amsmath"
    "CJK")
   (TeX-add-symbols
    '("referee" 1)
    '("YG" 1)
    '("todo" 1)
    '("Secref" 1)
    '("Tabref" 1)
    '("Figref" 1)
    '("Eqref" 1)
    '("code" 1)
    "mesa"
    "MESA"
    "kms"
    "kev"
    "gk"
    "zoph"
    "Msun"
    "Lsun")
   (LaTeX-add-labels
    "sec:intro"
    "sec:methods"
    "sec:best_model"
    "fig:HRD_both"
    "sec:MT"
    "fig:MT"
    "sec:rot"
    "sec:surf_rot"
    "fig:surf_rot"
    "sec:rot_comparison_single"
    "fig:struct_rot"
    "sec:mixing"
    "fig:D_mix"
    "fig:rho"
    "fig:n14"
    "sec:mix_comparison_single"
    "sec:surf_comp"
    "sec:discussion"
    "sec:single_star_uncertainties"
    "sec:bin_param"
    "sec:SN_comp"
    "sec:conclusions"
    "sec:software"
    "sec:res_tests"
    "fig:sp_test"
    "sec:X_fig"
    "fig:composition_huge"
    "fig:composition_movie"
    "fig:mix_movie")
   (LaTeX-add-bibliographies
    "./zeta_ophiuchi.bib"))
 :latex)

