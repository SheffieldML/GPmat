#! /bin/sh

if [ "A$1" = "A" ] ; then
    Rbinary="R" ;
else
    Rbinary="$1" ;
fi

ln -s tiger.Rnw.ignore tiger.Rnw
$Rbinary CMD Sweave tiger.Rnw
rm tiger.Rnw
mv tiger.tex tiger.Rnw
