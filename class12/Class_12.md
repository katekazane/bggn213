Class 12: Structural Bioinformatics and Drug Discovery
================

Cleaning up the Protein Target Structure
----------------------------------------

Start by downloading the target (the protein receptor) from the main PDB database. This time we'll be using the PDB ID "1hsg" the HIV protease.

Here we're writing a code that will allow us to switch out the protein (why we have the pdb.code defined).

``` r
library(bio3d)

pdb.code <- "1hsg"
file.name <- get.pdb(pdb.code)
```

    ## Warning in get.pdb(pdb.code): ./1hsg.pdb exists. Skipping download

This file is very large, so we want to extract the protein only segments. We'll then change it to a new PDB format file. We'll also be doing this for the bound ligand (will allow us to manipulate both).

``` r
hiv <- read.pdb(file.name)
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

Initially, we'll want just the protein. When you run it, check the non protein values (should be zero).

``` r
prot <- trim.pdb(hiv, "protein")
prot
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "protein")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1514,  XYZs#: 4542  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 0  (residues: 0)
    ##      Non-protein/nucleic resid values: [ none ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

We give the file a generic name to avoid hardocding 1hsg (keeps our code versatile).

``` r
prot.filename <- paste(pdb.code, "_protein.pdb", sep = "")
write.pdb(prot, file = prot.filename)
```

Doing the same for the ligand.

``` r
lig <- trim.pdb(hiv, "ligand")
lig
```

    ## 
    ##  Call:  trim.pdb(pdb = hiv, "ligand")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 45,  XYZs#: 135  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 45  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, helix, sheet, seqres, xyz,
    ##         calpha, call

``` r
lig.filename <- paste(pdb.code, "_ligand.pdb", sep = "")
write.pdb(lig, file = lig.filename)
```

Using The AutoDock Tools Package
--------------------------------

Once we open the 1hsg\_protein.pdb file, we add back the hydrogens and save the file as a .pdbqt Opening it we can see more of the electrostatic charges. Adding the hydrogens back gives us a lot of room to play around. It also assigns the charges to each molecule.

We can do the exact same thing for the ligands.

Converting Docking Results for Viewing in VMD
---------------------------------------------

``` r
res <- read.pdb("all.pdbqt", multi = TRUE)
res
```

    ## 
    ##  Call:  read.pdb(file = "all.pdbqt", multi = TRUE)
    ## 
    ##    Total Models#: 14
    ##      Total Atoms#: 50,  XYZs#: 2100  Chains#: 1  (values: B)
    ## 
    ##      Protein Atoms#: 0  (residues/Calpha atoms#: 0)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 50  (residues: 1)
    ##      Non-protein/nucleic resid values: [ MK1 (1) ]
    ## 
    ## + attr: atom, xyz, calpha, call

``` r
write.pdb(res, file = "results.pdb")
```

We can also compare docking results in R using RMSD (Root Mean Square Distance)

``` r
#We already ran the first part of this function.
#res <- read.pdb("all.pdbqt", multi=TRUE)
ori <- read.pdb("ligand.pdbqt")
rmsd(ori, res)
```

    ##  [1]  0.590 11.163 10.531  4.364 11.040  3.682  5.741  3.864  5.442 10.920
    ## [11]  4.318  6.249 11.084  8.929

This will compare all of the pairs (14x14 table).

``` r
rmsd(res)
```

    ## Warning in rmsd(res): No indices provided, using the 50 non NA positions

    ##         [,1]   [,2]   [,3]   [,4]   [,5]   [,6]   [,7]   [,8]   [,9]
    ##  [1,]  0.000 11.059 10.418  4.340 10.950  3.600  5.678  3.648  5.453
    ##  [2,] 11.059  0.000  4.168 10.675  4.924 10.902 11.246 10.473 10.804
    ##  [3,] 10.418  4.168  0.000 10.351  3.630 10.283 10.927  9.863 11.412
    ##  [4,]  4.340 10.675 10.351  0.000 10.435  5.544  4.239  5.763  6.970
    ##  [5,] 10.950  4.924  3.630 10.435  0.000 10.699 10.717 10.463 12.007
    ##  [6,]  3.600 10.902 10.283  5.544 10.699  0.000  4.827  2.867  6.335
    ##  [7,]  5.678 11.246 10.927  4.239 10.717  4.827  0.000  4.844  8.021
    ##  [8,]  3.648 10.473  9.863  5.763 10.463  2.867  4.844  0.000  6.296
    ##  [9,]  5.453 10.804 11.412  6.970 12.007  6.335  8.021  6.296  0.000
    ## [10,] 10.818  4.374  2.056 10.262  2.615 10.617 10.770 10.285 11.819
    ## [11,]  4.173 10.617 10.019  5.599 10.214  2.283  4.624  3.187  6.755
    ## [12,]  6.160 10.535 11.102  7.553 11.686  5.555  7.096  5.411  3.305
    ## [13,] 11.044  9.866 10.250 12.219 11.952 11.653 13.216 10.798  9.228
    ## [14,]  8.998 12.535 12.140  6.498 11.093  9.090  7.382  9.805 11.146
    ##        [,10]  [,11]  [,12]  [,13]  [,14]
    ##  [1,] 10.818  4.173  6.160 11.044  8.998
    ##  [2,]  4.374 10.617 10.535  9.866 12.535
    ##  [3,]  2.056 10.019 11.102 10.250 12.140
    ##  [4,] 10.262  5.599  7.553 12.219  6.498
    ##  [5,]  2.615 10.214 11.686 11.952 11.093
    ##  [6,] 10.617  2.283  5.555 11.653  9.090
    ##  [7,] 10.770  4.624  7.096 13.216  7.382
    ##  [8,] 10.285  3.187  5.411 10.798  9.805
    ##  [9,] 11.819  6.755  3.305  9.228 11.146
    ## [10,]  0.000 10.255 11.518 11.311 11.437
    ## [11,] 10.255  0.000  5.833 12.207  8.346
    ## [12,] 11.518  5.833  0.000  9.607 11.446
    ## [13,] 11.311 12.207  9.607  0.000 16.710
    ## [14,] 11.437  8.346 11.446 16.710  0.000

Ligand Based Approaches- Normal Mode Analysis
---------------------------------------------
