#! /bin/bash


for i in TRAV TRAJ TRAC TRBV TRBD TRBJ TRBC TRGV TRGJ TRGC TRDV TRDD TRDJ TRDC; do 
    wget --output-document=$i.web http://www.imgt.org/genedb/GENElect?query=7.2+$i\&species=Homo+sapiens
    #cat $i.web | perl -ne 'print if /^\>.+\|Homo sapiens\|.+$|^[ctag]+$/' | sed 's/([ctag])\n([ctag])//D' 
    cat $i.web | perl -ne 'print if /^\>.+\|Homo sapiens\|.+$|^[ctag]+$/' > $i.fasta
    rm -f $i.web
done
