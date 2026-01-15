function plink2_pca_approx(input_prefix; npcs=20, output_prefix=input_prefix)
    run(`plink2 \
            --bfile $input_prefix \
            --freq counts \
            --pca $npcs approx allele-wts vcols=chrom,ref,alt \
            --seed 1 \
            --out $output_prefix`
    )
    return CSV.read(string(output_prefix, ".eigenvec"), DataFrame)
end