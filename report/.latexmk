# Usage: latexmk report.tex
# Inspired by https://github.com/cdlm/infrastructure/blob/master/dotfiles/latexmkrc
$pdf_mode = 1;
$pdflatex = 'xelatex --shell-escape -file-line-error -halt-on-error -synctex=1 %O %S';

# http://tex.stackexchange.com/questions/1226/how-to-make-latexmk-use-makeglossaries
add_cus_dep('glo', 'gls', 0, 'run_makeglossaries');
add_cus_dep('acn', 'acr', 0, 'run_makeglossaries');
push @generated_exts, 'glo', 'gls', 'glg';
push @generated_exts, 'acn', 'acr', 'alg';
$clean_ext .= ' %R.ist %R.xdy';

sub run_makeglossaries {
  if ( $silent ) {
    system "makeglossaries -q '$_[0]'";
  }
  else {
    system "makeglossaries '$_[0]'";
  };
}
