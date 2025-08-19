export PATH=$PATH:$PWD/build
cd build
make -j
cd ../

MAT="../MAT_WEPP/SARS2-WBE/manuscript_art_unseen_10_dec_2023_dist_7/unseen_pruned_public-2023-12-25.all.masked.pb.gz"
sequences="../MAT_WEPP/SARS2-WBE/manuscript_art_unseen_10_dec_2023_dist_7/haplotype_removal.txt"
vcf="../MAT_WEPP/SARS2-WBE/manuscript_art_unseen_10_dec_2023_dist_7/my_vcf_samples.vcf"

ripples_server -i ${MAT} -s ${sequences} -v ${vcf} -T 32