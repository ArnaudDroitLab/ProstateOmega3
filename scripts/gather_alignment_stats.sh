grep -F 'processed' output/*/*.stderr |
    perl -ne '$_ =~ s/output\/(.*)_BP.* processed (.*) reads, (.*) reads pseudoaligned/$1\t$2\t$3/; print $_' |
    sed -e 's/,//g' > output/kallisto_stats.txt