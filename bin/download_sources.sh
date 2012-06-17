#!/bin/bash

# Downloads database sources from vendors

bin_dir="$( cd "$( dirname "$0" )" && pwd )"

if [ ! -d $bin_dir/../databases2 ] 
then
    mkdir -p $bin_dir/../databases2
fi

cd $bin_dir/../databases2

###
# Homo Sapiens
###

if [ ! -d ./hs ] 
then
    mkdir -p ./hs
fi

cd ./hs

### Localizations

if [ ! -d ./localization ] 
then
    mkdir -p ./localization
fi

cd ./localization

# Bacello
if [ ! -f pred_homo ]
then
    wget http://gpcr.biocomp.unibo.it/bacello/datasets/H.sapiens_annotations.zip -O ./bacello.zip
    unzip ./bacello.zip
    rm bacello.zip
fi   

# eSLDB
if [ ! -f eSLDB_Homo_sapiens.txt ]
then
    wget --post-data 'organism=Homo+sapiens&checkbox1_1=SI&checkbox1_2=SI&checkbox1_8=SI&f_type=plain&download=Download' \
    http://gpcr.biocomp.unibo.it/cgi-bin/predictors/esldb/full_download.cgi \
    -O ./eSLDB_Homo_sapiens.txt
fi

# GO
if [ ! -f mart_export.txt ]
then
    echo "HS GO database missing, please upload it manually"
fi    

# Human proteinpedia
if [ ! -f HUPA_Localization_Data.txt ]
then
    wget http://www.humanproteinpedia.org/HuPA_Download/FULL/HUPA_Localization_Data.txt.gz -O ./hupa.gz
    gzip -d -N ./hupa.gz
    rm ./hupa.gz
fi


# LOCATE
if [ ! -f LOCATE_human_v6_20081121.xml ]   
then
    wget http://locate.imb.uq.edu.au/info_files/LOCATE_human_v6_20081121.xml.zip -O ./locatexml.zip
    unzip locatexml.zip
    rm locatexml.zip
fi   

# MatrixDB Loc
if [ ! -f ECM_Protein_list_20100825.txt ]
then
    wget http://matrixdb.ibcp.fr/current/ECM_Protein_list_20100825.txt
fi

if [ ! -f Secreted_Protein_list_20100825.txt ]
then
    wget http://matrixdb.ibcp.fr/current/Secreted_Protein_list_20100825.txt
fi

if [ ! -f Membrane_Protein_list_20100825.txt ]
then
    wget http://matrixdb.ibcp.fr/current/Membrane_Protein_list_20100825.txt
fi

# Organelle
if [ ! -f organelle_hs.txt ]
then
    wget "http://organelledb.lsi.umich.edu/search-result.php?local=all&endoplasmic+reticulum=all&membrane+protein=all&miscellaneous+others=all&mitochondrion=all&nucleus=all&protein+complex=all&organism=human&searchtype=full&geneAdv=&flat_file=Download+in+Flat+File" -O ./organelle_hs.txt
fi  

# Pagosub
if [ ! -f homo_sapiens.csv ]
then
    wget http://webdocs.cs.ualberta.ca/%7Ebioinfo/PA/GOSUB/data/animal/homo_sapiens_csv.zip -O ./pagosub.zip
    unzip pagosub.zip
    rm pagosub.zip
fi

# Human Protein Atlas
if [ ! -f subcellular_location.csv ]
then
    wget http://www.proteinatlas.org/download/subcellular_location.csv.zip -O ./proteinAtlas.zip
    unzip proteinAtlas.zip
    rm proteinAtlas.zip
fi

### Interactions

if [ ! -d ../interaction ] 
then
    mkdir -p ../interaction
fi

cd ../interaction

# Biogrid
# @see at the end of this file

# CCSB
if [ ! -f HI2_2011.tsv ]
then
    wget http://interactome.dfci.harvard.edu/H_sapiens/download/HI2_2011.tsv
fi   

# ConsensusPathDB
if [ ! -f ConsensusPathDB_human_PPI ]
then
    wget http://cpdb.molgen.mpg.de/download/ConsensusPathDB_human_PPI.gz -O ./consensus.gz
    gzip -d -N ./consensus.gz
    rm ./consensus.gz
fi

# DiP
# @see at the end of this file

# HomoMINT
if [ ! -f homomint-full.mitab25.txt-binary.mitab26.txt ]
then
    wget ftp://mint.bio.uniroma2.it/pub/release/homomint/homomint-full.mitab25.txt-binary.mitab26.txt
fi

# HPRD
if [ ! -f BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt ]
then
    # download zip
    wget \
    -U "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.14) Gecko/20080418 Ubuntu/7.10 (gutsy) Firefox/2.0.0.14" \
    http://www.hprd.org/edownload/HPRD_Release9_041310 \
    -O ./hprd.tar.gz
    
    tar xvfz hprd.tar.gz HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt
    
    mv HPRD_Release9_062910/BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt ./BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt
    rm -r HPRD_Release9_062910
    
    # cleanup
    rm hprd.tar.gz
fi

# IntAct
# @see at the end of this file   

# Matrix
if [ ! -f MatrixDB_20120420.txt ]
then
    wget http://matrixdb.ibcp.fr/cgi-bin/MatrixDB_20120420.txt.gz -O ./matrix.gz
    gzip -d -N ./matrix.gz
    rm ./matrix.gz
fi

# MINT
if [ ! -f 2012-02-06-mint-human-binary.mitab26.txt ]
then
    wget "ftp://mint.bio.uniroma2.it/pub/release/mitab26/2012-02-06/2012-02-06-mint-human-binary.mitab26.txt"
fi

# MIPS
if [ ! -f allppis.xml ]
then
    wget http://mips.helmholtz-muenchen.de/proj/ppi/data/mppi.gz -O ./mips.gz
    gzip -d -N ./mips.gz
    rm ./mips.gz
fi

# PIP's
if [ ! -f PredictedInteractions1000.txt ]
then
    wget http://www.compbio.dundee.ac.uk/www-pips/downloads/PredictedInteractions1000.txt
fi 

cd ../../

###
# Drosophila Melanogaster 
###

if [ ! -d ./dm ] 
then
    mkdir -p ./dm
fi

cd ./dm

### Localizations

if [ ! -d ./localization ] 
then
    mkdir -p ./localization
fi

cd ./localization

# eSLDB
if [ ! -f Drosophila_melanogaster.BDGP5.4.49.pep.all.fa_subloc ]
then
    wget http://gpcr.biocomp.unibo.it/esldb/update/Drosophila_melanogaster.BDGP5.4.49.pep.all.fa_subloc
fi  

# GO
if [ ! -f mart_export.txt ]
then
    echo "DM GO database missing, please upload it manually"
fi  

# PaGOSUB
if [ ! -f drosophila_melanogaster.csv ]
then
    wget http://webdocs.cs.ualberta.ca/%7Ebioinfo/PA/GOSUB/data/animal/drosophila_melanogaster_csv.zip -O ./pagosub.zip
    unzip ./pagosub.zip
    rm ./pagosub.zip
fi   

# Organelle
if [ ! -f organelle_dm.txt ]
then
    wget "http://organelledb.lsi.umich.edu/search-result.php?local=all&endoplasmic+reticulum=all&membrane+protein=all&miscellaneous+others=all&mitochondrion=all&nucleus=all&protein+complex=all&organism=fly&searchtype=full&geneAdv=&flat_file=Download+in+Flat+File" -O ./organelle_dm.txt
fi 

### Interactions

if [ ! -d ../interaction ] 
then
    mkdir -p ../interaction
fi

cd ../interaction 

# Biogrid
# @see at the end of this file 

# DiP
# @see at the end of this file

# DroID
if [ ! -f finley_yth.txt ]
then
    wget http://www.droidb.org/data/Droid_2012_04/finley_yth.txt
fi

if [ ! -f CURAGEN_YTH.txt ]
then
    wget http://www.droidb.org/data/Droid_2012_04/CURAGEN_YTH.txt
fi

if [ ! -f HYBRIGENICS_YTH.txt ]
then
    wget http://www.droidb.org/data/Droid_2012_04/HYBRIGENICS_YTH.txt
fi

if [ ! -f DPIM_COAPCOMPLEX.txt ]
then
    wget http://www.droidb.org/data/Droid_2012_04/DPIM_COAPCOMPLEX.txt
fi

if [ ! -f FLY_OTHER_PHYSICAL.txt ]
then
    wget http://www.droidb.org/data/Droid_2012_04/FLY_OTHER_PHYSICAL.txt
fi

# IntAct
# @see at the end of this file

# MINT
if [ ! -f 2012-02-06-mint-Drosophila-binary.mitab26.txt ]
then
    wget "ftp://mint.bio.uniroma2.it/pub/release/mitab26/2012-02-06/2012-02-06-mint-Drosophila-binary.mitab26.txt"
fi

cd ../../

###
# Caenorhabditis elegans
###

if [ ! -d ./ce ] 
then
    mkdir -p ./ce
fi

cd ./ce

### Localizations

if [ ! -d ./localization ] 
then
    mkdir -p ./localization
fi

cd ./localization

# Bacello
if [ ! -f pred_cel ]
then
    wget http://gpcr.biocomp.unibo.it/bacello/datasets/C.elegans_annotations.zip -O ./bacello.zip
    unzip ./bacello.zip
    rm bacello.zip
fi

# eSLDB
if [ ! -f eSLDB_Caenorhabditis_elegans.txt ]
then
    wget --post-data 'organism=Caenorhabditis+elegans&checkbox1_1=SI&checkbox1_2=SI&checkbox1_8=SI&f_type=plain&download=Download' \
    http://gpcr.biocomp.unibo.it/cgi-bin/predictors/esldb/full_download.cgi \
    -O ./eSLDB_Caenorhabditis_elegans.txt
fi

# GO
if [ ! -f mart_export.txt ]
then
    echo "CE GO database missing, please upload it manually"
fi

# Organelle
if [ ! -f organelle_ce.txt ]
then
    wget "http://organelledb.lsi.umich.edu/search-result.php?local=all&endoplasmic+reticulum=all&membrane+protein=all&miscellaneous+others=all&mitochondrion=all&nucleus=all&protein+complex=all&organism=worm&searchtype=full&geneAdv=&flat_file=Download+in+Flat+File" -O ./organelle_ce.txt
fi

# Pagosub
if [ ! -f caenorhabditis_elegans.csv ]
then
    wget http://webdocs.cs.ualberta.ca/%7Ebioinfo/PA/GOSUB/data/animal/caenorhabditis_elegans_csv.zip -O ./pagosub.zip
    unzip pagosub.zip
    rm pagosub.zip
fi

### Interactions

if [ ! -d ../interaction ] 
then
    mkdir -p ../interaction
fi

cd ../interaction

# Biogrid
# @see at the end of this file

# CCSB
if [ ! -f wi8.txt ]
then
    wget http://interactome.dfci.harvard.edu/C_elegans/graphs/sequence_edges/wi8.txt
fi

# DiP
# @see at the end of this file

# IntAct
# @see at the end of this file

# MINT
if [ ! -f 2012-02-06-mint-Caenorhabditis-binary.mitab26.txt ]
then
    wget "ftp://mint.bio.uniroma2.it/pub/release/mitab26/2012-02-06/2012-02-06-mint-Caenorhabditis-binary.mitab26.txt"
fi

cd ../../

###
# Saccaromicies Cervisae
###

if [ ! -d ./sc ] 
then
    mkdir -p ./sc
fi

cd ./sc

### Localizations

if [ ! -d ./localization ] 
then
    mkdir -p ./localization
fi

cd ./localization

# Bacello
if [ ! -f pred_sce ]
then
    wget http://gpcr.biocomp.unibo.it/bacello/datasets/S.cerevisiae_annotations.zip -O ./bacello.zip
    unzip ./bacello.zip
    rm bacello.zip
fi

# eSLDB
if [ ! -f eSLDB_Saccharomyces_cerevisiae.txt ]
then
    wget --post-data 'organism=Saccharomyces+cerevisiae&checkbox1_1=SI&checkbox1_2=SI&checkbox1_8=SI&f_type=plain&download=Download' \
    http://gpcr.biocomp.unibo.it/cgi-bin/predictors/esldb/full_download.cgi \
    -O ./eSLDB_Saccharomyces_cerevisiae.txt
fi

# GO
if [ ! -f mart_export.txt ]
then
    echo "SC GO database missing, please upload it manually"
fi  

# Organelle
if [ ! -f organelle_sc.txt ]
then
    wget "http://organelledb.lsi.umich.edu/search-result.php?local=all&endoplasmic+reticulum=all&membrane+protein=all&miscellaneous+others=all&mitochondrion=all&nucleus=all&protein+complex=all&organism=yeast&searchtype=full&geneAdv=&flat_file=Download+in+Flat+File" -O ./organelle_sc.txt
fi

# Pagosub
if [ ! -f saccharomyces_cerevisiae.csv ]
then
    wget http://webdocs.cs.ualberta.ca/%7Ebioinfo/PA/GOSUB/data/fungi/saccharomyces_cerevisiae_csv.zip -O ./pagosub.zip
    unzip pagosub.zip
    rm pagosub.zip
fi

# pTARGET
# TODO source site is down or moved

### Interactions

if [ ! -d ../interaction ] 
then
    mkdir -p ../interaction
fi

cd ../interaction

# Biogrid
# @see at the end of this file

# CCSB
if [ ! -f Y2H_union.txt ]
then
    wget http://interactome.dfci.harvard.edu/S_cerevisiae/download/Y2H_union.txt
fi

# ConsensusPathDB
if [ ! -f ConsensusPathDB_yeast_PPI ]
then
    wget http://cpdb.molgen.mpg.de/download/ConsensusPathDB_yeast_PPI.gz -O ./consensus.gz
    gzip -d -N ./consensus.gz
    rm ./consensus.gz
fi

# DiP
# @see at the end of this file

# IntAct
# @see at the end of this file

# MINT
if [ ! -f 2012-02-06-mint-Saccharomyces-binary.mitab26.txt ]
then
    wget "ftp://mint.bio.uniroma2.it/pub/release/mitab26/2012-02-06/2012-02-06-mint-Saccharomyces-binary.mitab26.txt"
fi

cd ../../

# Biogrid
if [ ! -f ./hs/interaction/BIOGRID-ORGANISM-Homo_sapiens-3.1.81.tab2.txt ] \
|| [ ! -f ./dm/interaction/BIOGRID-ORGANISM-Drosophila_melanogaster-3.1.81.tab2.txt ] \
|| [ ! -f ./sc/interaction/BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.1.81.tab2.txt ] \
|| [ ! -f ./ce/interaction/BIOGRID-ORGANISM-Caenorhabditis_elegans-3.1.81.tab2.txt ]
then
    # download zip
    if [ ! -f BIOGRID-ORGANISM-3.1.81.tab2.zip ]
    then
        wget "http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.1.81/BIOGRID-ORGANISM-3.1.81.tab2.zip"
    fi
    
    # extract files
    unzip BIOGRID-ORGANISM-3.1.81.tab2.zip BIOGRID-ORGANISM-Homo_sapiens-3.1.81.tab2.txt BIOGRID-ORGANISM-Drosophila_melanogaster-3.1.81.tab2.txt BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.1.81.tab2.txt BIOGRID-ORGANISM-Caenorhabditis_elegans-3.1.81.tab2.txt
    
    mv BIOGRID-ORGANISM-Homo_sapiens-3.1.81.tab2.txt ./hs/interaction/
    mv BIOGRID-ORGANISM-Drosophila_melanogaster-3.1.81.tab2.txt ./dm/interaction/
    mv BIOGRID-ORGANISM-Saccharomyces_cerevisiae-3.1.81.tab2.txt ./sc/interaction/
    mv BIOGRID-ORGANISM-Caenorhabditis_elegans-3.1.81.tab2.txt ./ce/interaction/
    
    # remove zip
    rm BIOGRID-ORGANISM-3.1.81.tab2.zip     
fi

# DiP
if [ ! -f ./hs/interaction/Hsapi20120228.txt ] \
|| [ ! -f ./dm/interaction/Dmela20120228.txt ] \
|| [ ! -f ./sc/interaction/Scere20120228.txt ] \
|| [ ! -f ./ce/interaction/Celeg20120228.txt ]
then
    # login
    wget --save-cookies dip-cookies.txt \
     --keep-session-cookies \
     --post-data 'login=vAirDa99&pass=522QQ&lgn=1&Login=Login' \
     -U "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.14) Gecko/20080418 Ubuntu/7.10 (gutsy) Firefox/2.0.0.14" \
     http://dip.doe-mbi.ucla.edu/dip/Login.cgi \
     -O /dev/null
     
    # hs
    if [ ! -f ./hs/interaction/Hsapi20120228.txt ]
    then
        wget --load-cookies dip-cookies.txt \
        -U "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.14) Gecko/20080418 Ubuntu/7.10 (gutsy) Firefox/2.0.0.14" \
        http://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2012/tab25/Hsapi20120228.txt \
        -O ./hs/interaction/Hsapi20120228.txt
    fi
    
    # dm
    if [ ! -f ./dm/interaction/Dmela20120228.txt ]
    then
        wget --load-cookies dip-cookies.txt \
        -U "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.14) Gecko/20080418 Ubuntu/7.10 (gutsy) Firefox/2.0.0.14" \
        http://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2012/tab25/Dmela20120228.txt \
        -O ./dm/interaction/Dmela20120228.txt
    fi
    
    # sc
    if [ ! -f ./sc/interaction/Scere20120228.txt ]
    then
        wget --load-cookies dip-cookies.txt \
        -U "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.14) Gecko/20080418 Ubuntu/7.10 (gutsy) Firefox/2.0.0.14" \
        http://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2012/tab25/Scere20120228.txt \
        -O ./sc/interaction/Scere20120228.txt
    fi
    
    # ce
    if [ ! -f ./ce/interaction/Celeg20120228.txt ]
    then
        wget --load-cookies dip-cookies.txt \
        -U "Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.14) Gecko/20080418 Ubuntu/7.10 (gutsy) Firefox/2.0.0.14" \
        http://dip.doe-mbi.ucla.edu/dip/File.cgi?FN=2012/tab25/Celeg20120228.txt \
        -O ./ce/interaction/Celeg20120228.txt
    fi
    
    # cleanup
    rm dip-cookies.txt
fi

# IntAct
if [ ! -f ./hs/interaction/intact-hs.txt ] \
|| [ ! -f ./dm/interaction/intact-dm.txt ] \
|| [ ! -f ./sc/interaction/intact-sc.txt ] \
|| [ ! -f ./ce/interaction/intact-ce.txt ]
then
    if [ ! -f ./intact.zip ]
    then
        wget ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psimitab/intact.zip -O ./intact.zip
        unzip intact.zip
    fi     
    
    grep "taxid:9606(Human)	taxid:9606(Human)" intact.txt > ./hs/interaction/intact-hs.txt
    grep "taxid:7227(Fruit fly)	taxid:7227(Fruit fly)" intact.txt > ./dm/interaction/intact-dm.txt  
    grep "taxid:559292(Saccharomyces cerevisiae S288c)	taxid:559292(Saccharomyces cerevisiae S288c)" intact.txt > ./sc/interaction/intact-sc.txt
    grep "taxid:6239(Caenorhabditis elegans)	taxid:6239(Caenorhabditis elegans)" intact.txt > ./ce/interaction/intact-ce.txt
    
    # cleanup
    rm intact.txt
    rm intact.zip 
fi