#!/usr/bin/perl -w 
#
# Pad the first 6 columns

while (<>)
{
    # strip comments
    s/^([\*C].*$)//;
    # check lines with plain whitespace preceding a character on column 7
    s/^(\s{1,6})([a-zA-Z])/sprintf("%6s%1s",$1,$2)/e; 
    # check linrd with whitespace and numbers preceding a character on column 7
    s/^\s{0,3}([0-9]{1,5})\s{0,3}/sprintf("% 5d ",$1)/e;
    print $_;
}

