# Transposed VCF
Transposed VCF supplement MAT with information about ambiguous bases at sample nodes to aid tree optimization. All variants of a sample is written contagiously (in rows instead of columns), so new samples can be concatenated.
## Converting VCF to Transposed VCF
```
Transpose_vcf -v <input vcf> -o <output> -T <number of threads>
```
Output will be concatenated if exists. Currently, it can only utilize less than 5 threads, as it is throttled by decompression.
## Converting Transposed VCF to VCF
```
Transposed_vcf_to_vcf -i <transposed vcf> -o <vcf output> -T <number of threads>
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
        Compressed Data
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
                            <td> Position of mutation 1 (variant encoded) <\td>
                            <td> Position of mutation 2 (variant encoded) <\td>
                            <td> {allele at mutation 2 (one hot , 4 bits),allele at mutation 1 (one hot , 4 bits)} <\td>
                        </tr>
                        <tr> <td> ...... </td></tr>
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
                            <td> End of range (inclusive, variant encoded) <\td>
                            <td> Strat of range (inclusive, variant encoded, omitted if the same as end of range)  <\td>
                        </tr>
                        </tbody>
                        </table>
                    </td>
                </tr>
            </tbody>
        </table>
    </tr>
  </tbody>
</table>