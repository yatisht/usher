# Transposed VCF
Transposed VCF supplement MAT with information about ambiguous bases at sample nodes to aid tree optimization. All variants of a sample is written contagiously (in rows instead of columns), so new samples can be concatenated.
## Converting VCF to Transposed VCF
```
transpose_vcf -v <input vcf> -o <output> -T <number of threads>
```
Output will be concatenated if exists. Currently, it can only utilize less than 5 threads, as it is throttled by decompression.
## Converting Transposed VCF to VCF
```
transposed_vcf_to_vcf -i <transposed vcf> -o <vcf output> -r <reference fasta file> -T <number of threads>
```
## Use transposed VCF for optimization
```
matOptimize -i <usher protobuf> -V <transposed VCF> -o <output usher protobuf> -r <radius (4-6) > -T <number of threads>
```
## File Format
<table>
  <tbody>
    <tr>
        <td>Length of compressed Block (8 bytes)</td>
    </tr>
    <tr>
        <td>
        Compressed Data (zlib)
        <table>
            <tbody>
                <tr>
                    <td> Sample name (null terminated) </td>
                </tr>
                <tr>
                    <td> Called mutations (null terminated)
                        <table>
                        <tbody>
                            <tr> 
                              <td> Position of mutation 1 (variant encoded) 
                            <tr> 
                              <td> Position of mutation 2 (variant encoded, omitted if not present) 
                            <tr> 
                              <td> {allele at mutation 2 (one hot , 4 bits),allele at mutation 1 (one hot , 4 bits)} 
                        </tbody>
                        </table>
                    </td>
                </tr>
                <tr>
                    <td>
                        Completely ambiguous mutations (Ns) (Null terminated)
                        <table>
                        <tbody>
                        <tr>
                            <td> End of range (inclusive, variant encoded) 
                        <tr>
                            <td> Strat of range (inclusive, variant encoded, omitted if the same as end of range)
                        </tbody>
                        </table>
                    </td>
                </tr>
            </tbody>
        </table>
    </tr>
  </tbody>
</table>
Contig name is not recorded, as coronavirus only have one contig. If we need to handle multiple contig, concatnate all contig as a long reference, and let position be the position in that concatnated long contig. No change to file format is necessary, just transpose_vcf will need to take faidx to map contig+position to posiition in the concatnated long reference.
