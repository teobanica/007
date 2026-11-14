#############################################################################
##
##  dempsey.g  --  the metabelian groups Gamma_{K,q} (Flat linearity of uniform groups)
##
##  Usage:  in Terminal, go to the folder containing this file, then type
##
##      gap -q -o 2g dempsey.g
##
##  Output, written as it goes (so nothing is lost if you stop it):
##      dempsey.log       one line per uniform structure found, raw data
##      dempsey.progress  one line per group done, with timings
##  The script can be stopped at any time with Ctrl-C (twice), then quit;
##
##  What it does: for each pair [K, q], it builds the group
##  Gamma_{K,q} = (Z_K * Z_K)_met / qA, of order K^2 q^((K-1)^2), as a pc group,
##  and records, as in the census scripts, the character degrees, the spectral
##  data of the generators, the kernels and class sizes, and whether the S_2
##  symmetry holds and the diagonal kernel is abelian (field 3 of the log
##  line is q).  The pairs are in increasing order of size; the last ones,
##  [4,3] (order 314928) and [5,2] (order 1638400), are the heaviest.
##
#############################################################################

UG_TASKS := [ [3,3], [4,2], [3,7], [3,11], [3,13], [4,3], [5,2] ];   # pairs [K,q]
UG_NAME := "dempsey";
UG_OUT := Concatenation(UG_NAME, ".log");
UG_PROG := Concatenation(UG_NAME, ".progress");

UGWriteLine := function(file, str)
  local s;
  s := OutputTextFile(file, true);
  SetPrintFormattingStatus(s, false);
  AppendTo(s, str, "\n");
  CloseStream(s);
end;

UGIsSymmetric := function(G, gens)
  local M, perms, p, hom;
  M := Length(gens);
  if M = 1 then return true; fi;
  perms := [ Concatenation([2, 1], [3..M]), Concatenation([2..M], [1]) ];
  for p in perms do
    hom := GroupHomomorphismByImages(G, G, gens, gens{p});
    if hom = fail then return false; fi;
    if not IsBijective(hom) then return false; fi;
  od;
  return true;
end;

UGDiagonalKernel := function(G, gens)
  local M, L;
  M := Length(gens);
  L := Concatenation(GeneratorsOfGroup(DerivedSubgroup(G)),
                     List([2..M], i -> gens[i] * gens[1]^-1));
  if Length(L) = 0 then return TrivialSubgroup(G); fi;
  return Subgroup(G, L);
end;

UGLine := function(n, i, e, G, gens, K, M)
  local tbl, irr, sizes, w, dims, eig, row, chi, g, j, m, s, kers,
        lam, sym, L;
  tbl := CharacterTable(G);
  irr := Irr(tbl);
  sizes := SizesConjugacyClasses(tbl);
  w := E(K);
  dims := List(irr, chi -> chi[1]);
  eig := [];
  for chi in irr do
    row := [];
    for g in gens do
      for j in [0..K-1] do
        s := Sum([0..K-1], m -> (g^m)^chi * w^(-j*m)) / K;
        if not IsInt(s) then
          Error("non-integer eigenvalue multiplicity");
        fi;
        Add(row, s);
      od;
    od;
    Add(eig, row);
  od;
  kers := List(irr, ClassPositionsOfKernel);
  L := UGDiagonalKernel(G, gens);
  if Size(L) * K <> Size(G) then
    Error("diagonal kernel has wrong index");
  fi;
  if IsAbelian(L) then lam := 1; else lam := 0; fi;
  if UGIsSymmetric(G, gens) then sym := 1; else sym := 0; fi;
  return Concatenation("U|", String(n), "|", String(i), "|", String(e), "|",
    String(K), "|", String(M), "|", String(sym), "|", String(lam), "|",
    String(dims), "|", String(eig), "|", String(kers), "|", String(sizes));
end;


##  Gamma_{K,q} = (Z_K * Z_K)_met / qA, given by generators g1, g2 and relations
##  g1^K = g2^K = 1, c^q = 1 and [c,c'] = 1 for the commutators c = [g1^a, g2^b],
##  which form a free basis of the derived subgroup of Z_K * Z_K.

UGMetab := function(task)
  local K, q, d, n, F, g1, g2, cs, rels, a, b, P, epi, G, gens, t0;
  K := task[1]; q := task[2]; d := (K-1)^2; n := K^2 * q^d;
  t0 := Runtime();
  F := FreeGroup(2); g1 := F.1; g2 := F.2;
  cs := [];
  for a in [1..K-1] do for b in [1..K-1] do Add(cs, Comm(g1^a, g2^b)); od; od;
  rels := [g1^K, g2^K];
  Append(rels, List(cs, c -> c^q));
  for a in [1..Length(cs)] do for b in [a+1..Length(cs)] do
    Add(rels, Comm(cs[a], cs[b]));
  od; od;
  P := F / rels;
  epi := EpimorphismSolvableQuotient(P, n);
  G := Image(epi);
  if Size(G) <> n then
    UGWriteLine(UG_PROG, Concatenation("K=", String(K), " q=", String(q),
      " WRONG ORDER ", String(Size(G)), " instead of ", String(n)));
    return;
  fi;
  gens := List(GeneratorsOfGroup(P), x -> Image(epi, x));
  UGWriteLine(UG_OUT, UGLine(n, q, 0, G, gens, K, 2));
  UGWriteLine(UG_PROG, Concatenation("K=", String(K), " q=", String(q),
    " order=", String(n), " done, time=", String(Runtime() - t0), "ms"));
  Print("K=", K, " q=", q, " : done\n");
end;

UGWriteLine(UG_PROG, Concatenation("start ", UG_NAME));
for UG_TASK in UG_TASKS do
  UGMetab(UG_TASK);
od;
UGWriteLine(UG_PROG, "finished");
Print("finished\n");
QUIT;
