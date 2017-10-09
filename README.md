# Bioperl examples


# gfxx.pl sequence of scripts

After some time checking various gff utilities, what I found was that there wasn't a useful "getting to know your annotation" type script, one that would layout the information in a digestable and navigable manner. For that level of customisation, I endded up with the big platforms biopython and bioperl. But, and this has happened a few times, the former lags behind the latter, so I ended up witht he grandaddy bioperl.


## gf0.pl
taken from a web snippet.

## gf4.pl

This is quite a mature script now, and is able to deal with the S288 quite well. For the W303 ral, some work stil is needed .. the gene names seem to all have _W303 appended. Of course the RAL has very many contigs so will be trickier.


# Notes
Quotations not always necessary on value pair syntax
i.e. -asequence => $mysq will do fine. -asequence does not need to be in quotations.

You are better off not reducing the seqs to strings, and sicking to bioperl seq objects.

