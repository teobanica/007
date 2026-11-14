#############################################################################
##
##  cheetah.g  --  census of finite uniform groups (Flat linearity of uniform groups)
##
##  Usage:  in Terminal, go to the folder containing this file, then type
##
##      gap -q -o 3g cheetah.g
##
##  Output, written as it goes (so nothing is lost if you stop it):
##      cheetah.log       one line per uniform structure found, raw data
##      cheetah.progress  one line per order done, with timings
##  The script can be stopped at any time with Ctrl-C (twice), then quit;
##
##  What it does: for each task [K, M, nmin, nmax], it runs over the groups
##  G = SmallGroup(n, i) with nmin <= n <= nmax and K^M dividing n, keeps
##  those with abelianization Z_K^M, and finds all epimorphisms
##  Z_K^{*M} -> G up to automorphisms (GQuotients).  For each one, with
##  generators g_1..g_M, it records:
##    sym   : whether the permutations of g_1..g_M extend to automorphisms
##    lam   : whether the diagonal kernel ker(G -> Z_K, g_i -> 1) is abelian
##    dims  : the irreducible character degrees
##    eig   : for each irreducible, the multiplicities of E(K)^j, j=0..K-1,
##            as eigenvalues of g_1, ..., g_M (M*K integers per irreducible)
##    kers  : for each irreducible, the classes in its kernel
##    sizes : the class sizes
##  The flat covers, flat linearity and stationarity are then computed from
##  this data by census.py.
##
#############################################################################

UG_TASKS := [ [6,2,721,2000], [9,2,81,2000], [8,2,64,2000], [4,2,1001,2000], [3,3,730,2000] ];
UG_NAME := "cheetah";
UG_OUT := Concatenation(UG_NAME, ".log");
UG_PROG := Concatenation(UG_NAME, ".progress");
UG_MAXNR := 2500;     # skip orders having more groups than this

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

UGRunTask := function(task)
  local K, M, nmin, nmax, target, F, Fp, n, nr, i, G, quos, e, epi, gens,
        found, t0;
  K := task[1]; M := task[2]; nmin := task[3]; nmax := task[4];
  target := AbelianInvariants(DirectProduct(List([1..M], i -> CyclicGroup(K))));
  F := FreeGroup(M);
  Fp := F / List(GeneratorsOfGroup(F), x -> x^K);
  for n in [nmin..nmax] do
    if n mod K^M = 0 then
      nr := NrSmallGroups(n);
      if nr > UG_MAXNR then
        UGWriteLine(UG_OUT, Concatenation("SKIP|", String(n), "|", String(nr),
                                     "|", String(K), "|", String(M)));
        UGWriteLine(UG_PROG, Concatenation("K=", String(K), " M=", String(M),
                  " n=", String(n), " skipped (", String(nr), " groups)"));
      else
        t0 := Runtime();
        found := 0;
        for i in [1..nr] do
          G := SmallGroup(n, i);
          if AbelianInvariants(G) = target then
            quos := GQuotients(Fp, G);
            e := 0;
            for epi in quos do
              e := e + 1;
              gens := List(GeneratorsOfGroup(Fp), x -> Image(epi, x));
              if ForAll(gens, g -> Order(g) = K) then
                UGWriteLine(UG_OUT, UGLine(n, i, e, G, gens, K, M));
                found := found + 1;
              fi;
            od;
          fi;
        od;
        UGWriteLine(UG_PROG, Concatenation("K=", String(K), " M=", String(M),
                  " n=", String(n), " groups=", String(nr),
                  " uniform structures=", String(found),
                  " time=", String(Runtime() - t0), "ms"));
        Print("K=", K, " M=", M, " n=", n, " : ", found,
              " uniform structures\n");
      fi;
    fi;
  od;
end;

UGWriteLine(UG_PROG, Concatenation("start ", UG_NAME));
for UG_TASK in UG_TASKS do
  UGRunTask(UG_TASK);
od;
UGWriteLine(UG_PROG, "finished");
Print("finished\n");
QUIT;
