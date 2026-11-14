Ancillary files for "Flat linearity of uniform groups"

Python programs (Python 3, with numpy and scipy):

  flatlin.py      core library: flat representations, flat covers, flat
                  linearity, stationarity and dimension-type linear programs
  examples.py     basic families of uniform groups (section 4), output in
                  examples.log
  metab.py        the groups Gamma_{K,q}, from the normal form c^v g1^x g2^y
                  (section 4), run as "python3 metab.py 3,2 3,3", output in
                  metab.log
  bloch.py        flat representations of the reflection groups
                  Gamma_N^{[s]}, via the Bloch sphere (section 7), output in
                  bloch.log
  flatpart.py     the flat part of Gamma_N^{[4]} (section 7), output in
                  flatpart.log
  census.py       analysis of the GAP census logs, output in census/census.txt
  census_sym.py   faster variant of census.py, restricted to the symmetric
                  epimorphisms, used for liston.log and canelo.log, output in
                  census/liston.json and census/canelo.json
  k9.py           the uniform structures at K=9, M=2 from Proposition 8.6,
                  run as "python3 k9.py census/liston.log"
  explore.py      random search for uniform groups among permutation groups,
                  output in census/explore-*.log
  heights.py      flatness at height m (section 9), output in
                  census/heights.log

GAP programs (GAP 4.16.1, with the SmallGroups library):

  gap/template.g  template for the census programs
  gap/manul.g, gap/leopard.g, gap/ali.g, gap/cheetah.g, gap/frazier.g,
  gap/liston.g, gap/canelo.g
                  SmallGroups census (section 8), each covering a range of
                  K, M and group orders
  gap/dempsey.g   the groups Gamma_{K,q} of larger order, from their
                  presentation (section 4)

The raw output of each GAP program is in census/NAME.log, with its progress
record in census/NAME.progress. The largest log, census/ali.log.gz, is
gzipped, and census.py reads it as is.
