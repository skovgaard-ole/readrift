# Example data

A real dataset small enough to live in the repository: Oxford Nanopore (GridION)
reads of *Methylocystis iwaonis* SS37A-Re aligned to the strain's own genome.
It lets you try ReadRift straight after installing it, without BLAST, and it is
what `tests/test_example.py` checks ReadRift's results against.

| File | Size | Contents |
|---|---|---|
| `AP027148.gb.gz` | 2.8 MB | the reference: chromosome `AP027142`, plasmids `AP027143`–`AP027148`, and four short gene records of the same organism (16S rRNA, *pmoA*, *mmoX*, *mxaF*), with their annotations |
| `DRR325755_20x.btop.gz` | 23.4 MB | the alignment of the first 7,435 reads of SRA run `DRR325755`, in the format ReadRift reads, including the `btop` trace column |

## Try it

From the repository folder:

```bash
python -m readrift examples/AP027148.gb.gz examples/DRR325755_20x.btop.gz -x 20 --identity
```

It takes about a minute, most of it drawing the PDF. It writes
`DRR325755_20x.pdf` and `DRR325755_20x.readriftdb.npz` into the folder you run
it from, then explore the result with

```bash
python -m readrift browse DRR325755_20x.readriftdb.npz
```

You should see:

| | |
|---|---|
| Reads | 3,120 |
| Undivided | 2,663 (85.35%) |
| Divided, short distance | 340 (10.90%) |
| Divided, long distance | 117 (3.75%) |
| Reads indicating inversions | 12 |
| Mean coverage (matched bases) | 19.62× |
| Read N50 | 34,922 bp |
| Structural events | 318 |
| Fold-back reads dropped | 156 |
| Median alignment identity | 94.0% (with `--identity`) |

## Why 20×

The full run holds 817× coverage in a 985 MB alignment. `-x 20` stops reading
once 20× is reached, so a 20× run only ever reads the first 7,435 reads, and
those are exactly the reads included here. These numbers match a `-x 20` run on
the complete alignment read for read, not just in total.

## How it was made

With the commands the main README documents, on NCBI BLAST+ 2.13.0:

1. **Reference FASTA for the database.** `makeblastdb` reads only FASTA, so the
   11 sequences were written out of `AP027148.gb`, each named by its `VERSION`
   (`>AP027142.1` …). Downloading the same accessions as FASTA from NCBI gives
   the same sequences under the same names.
2. **Reads.** `DRR325755.1` to `DRR325755.7435` from the run's FASTA, in their
   original order.
3. **Alignment:**

   ```bash
   makeblastdb -in reference.fa -dbtype nucl -out refdb
   blastn -db refdb -query reads_20x.fa -out DRR325755_20x.btop -outfmt "6 qseqid sframe qstart qend sstart send qlen sseqid btop" -num_threads 7
   ```

   It took 4.5 minutes and wrote 448,764 lines. The line endings were
   converted to LF, as BLAST writes them on Linux and macOS, and the file was
   compressed with `gzip -9`.

The reads themselves are not included: even the 20× subset is 104 MB of
sequence. Fetch them from the Sequence Read Archive to repeat step 3.

## Credit and licence

The data are **not** covered by ReadRift's MIT licence. They are public
sequence data, deposited in the International Nucleotide Sequence Database
Collaboration (via DDBJ), which places no restrictions on their use or
redistribution. If you use them, cite the original work:

> Kaise H, Sawadogo JB, Alam MS, Ueno C, Dianou D, Shinjo R, Asakawa S (2023).
> *Methylocystis iwaonis* sp. nov., a type II methane-oxidizing bacterium from
> surface soil of a rice paddy field in Japan, and emended description of the
> genus *Methylocystis* (ex Whittenbury et al. 1970) Bowman et al. 1993.
> *Int J Syst Evol Microbiol* 73(6). doi:10.1099/ijsem.0.005925

BioProject PRJDB12481 · BioSample SAMD00412743 · SRA run DRR325755
