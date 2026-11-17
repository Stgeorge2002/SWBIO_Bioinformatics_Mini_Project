# Ubuntu bacterial-genome analysis commands

This runbook downloads short- and long-read data for the 70 *Klebsiella
pneumoniae* ST258 isolates selected from BioProject PRJNA1133668. It uses
individual Docker and command-line steps.

Use Ubuntu x86-64 with Docker, at least 12 CPU cores, 32 GB RAM, and enough disk
for the read data, databases, and results. Run the blocks in order from the
same Bash terminal.

## 1. Install prerequisites and create the work folder

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates curl wget pigz docker.io python3
sudo systemctl enable --now docker
sudo usermod -aG docker "$USER"
```

After logging back in:

```bash
docker info
export RUN_ROOT="$HOME/kp-st258-analysis"
mkdir -p "$RUN_ROOT"
cd "$RUN_ROOT"
```

Install SRA Toolkit 3.4.1:

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.4.1/sratoolkit.3.4.1-ubuntu64.tar.gz
tar -xzf sratoolkit.3.4.1-ubuntu64.tar.gz
export PATH="$RUN_ROOT/sratoolkit.3.4.1-ubuntu64/bin:$PATH"
fasterq-dump --version
```

## 2. Download short and long reads

Use the [NCBI SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1133668)
to select the 70 ST258 isolates and export the run accessions for their short-
and long-read libraries. Add one row per SRA run to this tab-separated file.
Use simple isolate IDs in the first column and `short` or `long` in the second.

```bash
export READS="$RUN_ROOT/reads/short"
export LONG_READS="$RUN_ROOT/reads/long"
export SRA_CACHE="$RUN_ROOT/sra-cache"
export SRA_TMP="$RUN_ROOT/sra-tmp"
mkdir -p "$READS" "$LONG_READS" "$SRA_CACHE" "$SRA_TMP"
printf 'sample\tread_type\taccession\n' > "$RUN_ROOT/sra_runs.tsv"
```

Fill `sra_runs.tsv` with the selected accessions, then download and convert them:

```bash
while IFS=$'\t' read -r sample read_type accession; do
  [ "$sample" = "sample" ] && continue
  [ -n "$sample" ] && [ -n "$read_type" ] && [ -n "$accession" ] || continue
  prefetch "$accession" --max-size u --output-directory "$SRA_CACHE"
  run_dir="$SRA_TMP/$accession"
  mkdir -p "$run_dir"
  if [ "$read_type" = "short" ]; then
    fasterq-dump "$SRA_CACHE/$accession/$accession.sra" \
      --split-files --threads 12 --outdir "$run_dir"
    pigz -f "$run_dir/${accession}_1.fastq" "$run_dir/${accession}_2.fastq"
    cat "$run_dir/${accession}_1.fastq.gz" >> "$READS/${sample}_R1.fastq.gz"
    cat "$run_dir/${accession}_2.fastq.gz" >> "$READS/${sample}_R2.fastq.gz"
  elif [ "$read_type" = "long" ]; then
    fasterq-dump "$SRA_CACHE/$accession/$accession.sra" \
      --threads 12 --outdir "$run_dir"
    pigz -f "$run_dir/${accession}.fastq"
    cat "$run_dir/${accession}.fastq.gz" >> "$LONG_READS/${sample}.fastq.gz"
  else
    echo "Unknown read type for $accession: $read_type" >&2
    exit 1
  fi
done < "$RUN_ROOT/sra_runs.tsv"

find "$READS" "$LONG_READS" -type f -name '*.fastq.gz' -print
```

The script combines multiple runs for the same isolate by appending compressed
FASTQ files. Keep the read folders empty before running it to avoid duplicate
data if the download is repeated.

## 3. Pull analysis images and initialize databases

```bash
for image in \
  staphb/fastp:1.4.0 \
  staphb/shovill:1.4.2 \
  staphb/unicycler:0.5.1 \
  staphb/bandage:0.9.0 \
  quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1 \
  staphb/bakta:1.12.2 \
  quay.io/biocontainers/panaroo:1.8.1--pyhdfd78af_0 \
  staphb/snippy:4.6.0 \
  ubuntu:24.04
do
  docker pull "$image"
done
```

Initialize the full Bakta database and MOB-suite database:

```bash
export BAKTA_DB="$RUN_ROOT/databases/bakta"
export MOB_DB="$RUN_ROOT/databases/mob_suite"
mkdir -p "$BAKTA_DB" "$MOB_DB"
docker run --rm --entrypoint /bin/bash -v "$BAKTA_DB:/db" staphb/bakta:1.12.2 \
  -lc "bakta_db download --output /db --type full"
docker run --rm -v "$MOB_DB:/db" quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1 \
  mob_init --database_directory /db
```

## 4. Build the x86 mge-cluster image

```bash
mkdir -p /tmp/mge-cluster-build
cat > /tmp/mge-cluster-build/Dockerfile <<'DOCKERFILE'
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl bzip2 procps \
 && rm -rf /var/lib/apt/lists/*
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C /tmp \
 && install -m 0755 /tmp/bin/micromamba /usr/local/bin/micromamba \
 && rm -rf /tmp/bin
RUN micromamba create -y -p /opt/env/mge \
    -c conda-forge -c bioconda --channel-priority strict \
    python=3.10 pip wheel 'unitig-caller=1.3.2' hdbscan opentsne pandas \
    'numpy>=1.16' 'scikit-learn>=1.0' joblib \
 && micromamba clean -y -a
RUN curl -Ls \
    https://gitlab.com/sirarredondo/mge-cluster/-/archive/1.1.0/mge-cluster-1.1.0.tar.gz \
    -o /tmp/mge-cluster.tar.gz \
 && mkdir -p /tmp/mge-cluster-src \
 && tar -xzf /tmp/mge-cluster.tar.gz -C /tmp/mge-cluster-src --strip-components=1 \
 && /opt/env/mge/bin/pip install --no-cache-dir 'setuptools==65.5.1' \
 && /opt/env/mge/bin/pip install --no-cache-dir --no-deps --no-build-isolation /tmp/mge-cluster-src \
 && rm -rf /tmp/mge-cluster-src /tmp/mge-cluster.tar.gz
ENV PATH="/opt/env/mge/bin:/usr/local/bin:/usr/bin:/bin"
ENV LD_LIBRARY_PATH="/opt/env/mge/lib"
RUN /opt/env/mge/bin/python -c 'import unitig_query, mge_cluster, hdbscan, openTSNE' \
 && test -x /opt/env/mge/bin/mge_cluster
DOCKERFILE

docker build -t mge-cluster:1.1.0 /tmp/mge-cluster-build
docker run --rm mge-cluster:1.1.0 mge_cluster --help >/dev/null
rm -rf /tmp/mge-cluster-build
```

## 5. Run short-read assembly, annotation, and plasmid detection

These analysis commands use short reads for Shovill assembly. Long reads are
downloaded and retained. To use both read types with hybrid assembly instead,
run Section 5B in place of this section. After either choice, continue with
Sections 6–8.

```bash
export RESULTS="$RUN_ROOT/cli-results"
export DOCKER_USER="$(id -u):$(id -g)"
mkdir -p "$RESULTS/samples" "$RESULTS/plasmids"

for r1 in "$READS"/*_R1.fastq.gz; do
  [ -f "$r1" ] || { echo "No *_R1.fastq.gz files in $READS" >&2; exit 1; }
  sample="$(basename "$r1" _R1.fastq.gz)"
  r2="$READS/${sample}_R2.fastq.gz"
  [ -s "$r2" ] || { echo "Missing mate: $r2" >&2; exit 1; }
  sample_dir="$RESULTS/samples/$sample"
  mkdir -p "$sample_dir"/{fastp,assembly,bandage,bakta,mob}

  # FASTP
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/fastp" \
    staphb/fastp:1.4.0 \
    fastp -i "$r1" -I "$r2" \
      -o "$sample_dir/fastp/${sample}_cleaned_R1.fastq.gz" \
      -O "$sample_dir/fastp/${sample}_cleaned_R2.fastq.gz" \
      -j "$sample_dir/fastp/${sample}_fastp.json" \
      -h "$sample_dir/fastp/${sample}_fastp.html" \
      -q 20 -u 20 -5 20 -W 4 \
      --cut_front --cut_tail \
      --cut_front_window_size 1 --cut_front_mean_quality 20 \
      --cut_tail_window_size 1 --cut_tail_mean_quality 20 \
      -l 50 -w 2

  # Shovill
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/assembly" \
    staphb/shovill:1.4.2 \
    shovill \
      --R1 "$sample_dir/fastp/${sample}_cleaned_R1.fastq.gz" \
      --R2 "$sample_dir/fastp/${sample}_cleaned_R2.fastq.gz" \
      --outdir "$sample_dir/assembly/shovill" \
      --force --cpus 12 --ram 27 --minlen 500 --mincov 5
  cp "$sample_dir/assembly/shovill/contigs.fa" "$sample_dir/assembly/${sample}.fasta"
  graph="$(find "$sample_dir/assembly/shovill" -type f \( -name '*.gfa' -o -name '*.fastg' \) -print -quit)"
  if [ -n "$graph" ]; then
    cp "$graph" "$sample_dir/assembly/${sample}.gfa"

    # Bandage
    docker run --rm --user "$DOCKER_USER" \
      -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/bandage" \
      staphb/bandage:0.9.0 \
      Bandage info "$sample_dir/assembly/${sample}.gfa" \
      > "$sample_dir/bandage/${sample}_bandage.info.txt"
    docker run --rm --user "$DOCKER_USER" \
      -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/bandage" \
      -e QT_QPA_PLATFORM=offscreen \
      staphb/bandage:0.9.0 \
      Bandage image "$sample_dir/assembly/${sample}.gfa" \
        "$sample_dir/bandage/${sample}_bandage.png" --height 1000 --width 1000
  else
    echo "No Shovill graph for $sample; Bandage skipped." >&2
  fi

  # MOB-recon
  docker run --rm --user root \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/mob" \
    quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1 bash -lc '
      db="'"$MOB_DB"'"
      egg_db="$(python3 -c '\''import mob_suite, os; print(os.path.join(os.path.dirname(mob_suite.__file__), "databases"))'\'')"
      mkdir -p "$egg_db"
      ln -sf "$db/taxa.sqlite" "$egg_db/taxa.sqlite"
      [ ! -s "$db/taxa.sqlite.traverse.pkl" ] || ln -sf "$db/taxa.sqlite.traverse.pkl" "$egg_db/taxa.sqlite.traverse.pkl"
      mob_recon --infile "'"$sample_dir"'/assembly/'"$sample"'.fasta" \
        --outdir "'"$sample_dir"'/mob/'"$sample"'_mob_recon" \
        --sample_id "'"$sample"'" --num_threads 8 --min_length 1000 \
        --database_directory "$db" --force
    '

  # Bakta on chromosome contigs
  chromosome="$sample_dir/mob/${sample}_mob_recon/chromosome.fasta"
  [ -s "$chromosome" ] || { echo "Missing MOB-suite chromosome FASTA: $chromosome" >&2; exit 1; }
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/bakta" \
    staphb/bakta:1.12.2 \
    bakta --output "$sample_dir/bakta/output" --prefix "$sample" \
      --db "$BAKTA_DB" --threads 12 --locus-tag STRPN \
      --min-contig-length 200 \
      --skip-plot --skip-trna --skip-tmrna --skip-rrna \
      --skip-ncrna --skip-crispr --force \
      "$chromosome"
done
```

## 5B. Hybrid assembly with Unicycler (alternative to Section 5)

Choose this section instead of Section 5, not in addition to it. It uses
quality-trimmed short-read pairs and each sample's long reads for Unicycler.
The existing Sections 6–8 then work with the hybrid assemblies. The separate
`RESULTS` directory keeps these results apart from short-read-only results.

```bash
export RESULTS="$RUN_ROOT/hybrid-results"
export DOCKER_USER="$(id -u):$(id -g)"
mkdir -p "$RESULTS/samples" "$RESULTS/plasmids"

for r1 in "$READS"/*_R1.fastq.gz; do
  [ -f "$r1" ] || { echo "No *_R1.fastq.gz files in $READS" >&2; exit 1; }
  sample="$(basename "$r1" _R1.fastq.gz)"
  r2="$READS/${sample}_R2.fastq.gz"
  long_reads="$LONG_READS/${sample}.fastq.gz"
  [ -s "$r2" ] || { echo "Missing short-read mate: $r2" >&2; exit 1; }
  [ -s "$long_reads" ] || { echo "Missing long reads: $long_reads" >&2; exit 1; }

  sample_dir="$RESULTS/samples/$sample"
  mkdir -p "$sample_dir"/{fastp,assembly,bandage,bakta,mob}

  # Quality-trim the paired short reads
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/fastp" \
    staphb/fastp:1.4.0 \
    fastp -i "$r1" -I "$r2" \
      -o "$sample_dir/fastp/${sample}_cleaned_R1.fastq.gz" \
      -O "$sample_dir/fastp/${sample}_cleaned_R2.fastq.gz" \
      -j "$sample_dir/fastp/${sample}_fastp.json" \
      -h "$sample_dir/fastp/${sample}_fastp.html" \
      -q 20 -u 20 -5 20 -W 4 \
      --cut_front --cut_tail \
      --cut_front_window_size 1 --cut_front_mean_quality 20 \
      --cut_tail_window_size 1 --cut_tail_mean_quality 20 \
      -l 50 -w 2

  # Hybrid assembly with trimmed paired reads and the long reads
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/assembly" \
    staphb/unicycler:0.5.1 \
    unicycler \
      -1 "$sample_dir/fastp/${sample}_cleaned_R1.fastq.gz" \
      -2 "$sample_dir/fastp/${sample}_cleaned_R2.fastq.gz" \
      -l "$long_reads" \
      -o "$sample_dir/assembly/unicycler" \
      --threads 12 --mode conservative --min_fasta_length 500
  cp "$sample_dir/assembly/unicycler/assembly.fasta" \
    "$sample_dir/assembly/${sample}.fasta"
  cp "$sample_dir/assembly/unicycler/assembly.gfa" \
    "$sample_dir/assembly/${sample}.gfa"

  # Bandage graph summary and image
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/bandage" \
    staphb/bandage:0.9.0 \
    Bandage info "$sample_dir/assembly/${sample}.gfa" \
    > "$sample_dir/bandage/${sample}_bandage.info.txt"
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/bandage" \
    -e QT_QPA_PLATFORM=offscreen \
    staphb/bandage:0.9.0 \
    Bandage image "$sample_dir/assembly/${sample}.gfa" \
      "$sample_dir/bandage/${sample}_bandage.png" --height 1000 --width 1000

  # Identify plasmids and annotate chromosome contigs
  docker run --rm --user root \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/mob" \
    quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1 bash -lc '
      db="'"$MOB_DB"'"
      egg_db="$(python3 -c '\''import mob_suite, os; print(os.path.join(os.path.dirname(mob_suite.__file__), "databases"))'\'')"
      mkdir -p "$egg_db"
      ln -sf "$db/taxa.sqlite" "$egg_db/taxa.sqlite"
      [ ! -s "$db/taxa.sqlite.traverse.pkl" ] || ln -sf "$db/taxa.sqlite.traverse.pkl" "$egg_db/taxa.sqlite.traverse.pkl"
      mob_recon --infile "'"$sample_dir"'/assembly/'"$sample"'.fasta" \
        --outdir "'"$sample_dir"'/mob/'"$sample"'_mob_recon" \
        --sample_id "'"$sample"'" --num_threads 8 --min_length 1000 \
        --database_directory "$db" --force
    '

  chromosome="$sample_dir/mob/${sample}_mob_recon/chromosome.fasta"
  [ -s "$chromosome" ] || { echo "Missing MOB-suite chromosome FASTA: $chromosome" >&2; exit 1; }
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$sample_dir/bakta" \
    staphb/bakta:1.12.2 \
    bakta --output "$sample_dir/bakta/output" --prefix "$sample" \
      --db "$BAKTA_DB" --threads 12 --locus-tag STRPN \
      --min-contig-length 200 \
      --skip-plot --skip-trna --skip-tmrna --skip-rrna \
      --skip-ncrna --skip-crispr --force \
      "$chromosome"
done
```

## 6. Annotate plasmids and run mge-cluster and Panaroo

```bash
mkdir -p "$RESULTS/plasmids/fastas" "$RESULTS/plasmids/gffs" \
  "$RESULTS/plasmids/families" "$RESULTS/panaroo" "$RESULTS/mge"

# Bakta on plasmids
while IFS= read -r fasta; do
  sample="$(basename "$(dirname "$(dirname "$(dirname "$fasta")")")")"
  plasmid="$(basename "$fasta")"
  stem="${sample}__${plasmid%.*}"
  cp "$fasta" "$RESULTS/plasmids/fastas/${stem}.fasta"
  mkdir -p "$RESULTS/plasmids/gffs/$stem"
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$RESULTS/plasmids/gffs/$stem" \
    staphb/bakta:1.12.2 \
    bakta --output "$RESULTS/plasmids/gffs/$stem" --prefix "$stem" \
      --db "$BAKTA_DB" --threads 4 --locus-tag STRPN \
      --min-contig-length 200 \
      --skip-plot --skip-trna --skip-tmrna --skip-rrna \
      --skip-ncrna --skip-crispr --force \
      "$RESULTS/plasmids/fastas/${stem}.fasta"
  cp "$RESULTS/plasmids/gffs/$stem/${stem}.gff3" \
    "$RESULTS/plasmids/gffs/${stem}.gff3"
done < <(find "$RESULTS/samples" -path '*/mob/*_mob_recon/plasmid_*.fasta' -type f | sort)

# Prepare mge-cluster inputs
mkdir -p "$RESULTS/mge/sequences"
: > "$RESULTS/mge/mge_input.txt"
for fasta in "$RESULTS"/plasmids/fastas/*.fasta; do
  [ -f "$fasta" ] || continue
  stem="$(basename "$fasta" .fasta)"
  awk -v id="$stem" 'BEGIN{n=0} /^>/{n++; print ">" id (n==1 ? "" : "_contig" n); next} {print}' \
    "$fasta" > "$RESULTS/mge/sequences/${stem}.fasta"
  printf '%s\n' "$RESULTS/mge/sequences/${stem}.fasta" >> "$RESULTS/mge/mge_input.txt"
done
plasmid_count="$(wc -l < "$RESULTS/mge/mge_input.txt")"
[ "$plasmid_count" -ge 2 ] || { echo "mge-cluster needs at least 2 plasmids" >&2; exit 1; }
min_cluster=4
perplexity=30
[ "$min_cluster" -le "$plasmid_count" ] || min_cluster="$plasmid_count"
[ "$perplexity" -lt "$plasmid_count" ] || perplexity=$((plasmid_count - 1))

# MGE-cluster
docker run --rm --user "$DOCKER_USER" \
  -v "$RUN_ROOT:$RUN_ROOT" -w "$RESULTS/mge" \
  mge-cluster:1.1.0 \
  mge_cluster --create --input "$RESULTS/mge/mge_input.txt" \
    --outdir "$RESULTS/mge/mge_cluster_out" --prefix plasmids \
    --kmer 31 --min_cluster "$min_cluster" --variance 0.01 \
    --perplexity "$perplexity" --threads 8
```

Convert the mge-cluster table into the family map, then assign plasmid GFFs:

```bash
python3 - "$RESULTS/mge/mge_cluster_out/plasmids_results.csv" \
  "$RESULTS/mge/mge_input.txt" "$RESULTS/mge/typing.tsv" <<'PY'
import csv, os, sys
results, manifest, output = sys.argv[1:]
expected = [os.path.splitext(os.path.basename(x.strip()))[0]
            for x in open(manifest) if x.strip()]
rows = {}
with open(results, newline="") as handle:
    for row in csv.DictReader(handle):
        rows[os.path.splitext(os.path.basename(row["Sample_Name"]))[0]] = row["Standard_Cluster"].strip()
with open(output, "w") as handle:
    handle.write("plasmid_id\tfamily\n")
    for stem in expected:
        cluster = rows[stem]
        token = cluster[:-2] if cluster.endswith(".0") else cluster
        family = "mge_na_" + stem.replace("__", "_") if token in {"-1", "-", "nan", "None", ""} else "mge_" + token
        handle.write(f"{stem}\t{family}\n")
PY

while IFS=$'\t' read -r plasmid family; do
  [ "$plasmid" = "plasmid_id" ] && continue
  mkdir -p "$RESULTS/plasmids/families/$family"
  cp "$RESULTS/plasmids/gffs/${plasmid}.gff3" "$RESULTS/plasmids/families/$family/"
done < "$RESULTS/mge/typing.tsv"

# Panaroo on chromosome GFFs
mapfile -t chromosome_gffs < <(find "$RESULTS/samples" -path '*/bakta/output/*.gff3' -type f | sort)
docker run --rm --user "$DOCKER_USER" \
  -v "$RUN_ROOT:$RUN_ROOT" -w "$RESULTS/panaroo" \
  quay.io/biocontainers/panaroo:1.8.1--pyhdfd78af_0 \
  panaroo -i "${chromosome_gffs[@]}" -o "$RESULTS/panaroo/chromosome" \
    --clean-mode sensitive --remove-invalid-genes --all_seq_in_graph \
    --len_dif_percent 0.95 --length_outlier_support_proportion 0.95 \
    --min_trailing_support 2 --trailing_recursive 1 \
    --min_edge_support_sv 0 --remove_by_consensus False \
    --edge_support_threshold 0 --threads 12

# Panaroo once per plasmid family
for family_dir in "$RESULTS"/plasmids/families/*; do
  [ -d "$family_dir" ] || continue
  mapfile -t family_gffs < <(find "$family_dir" -maxdepth 1 -name '*.gff3' -type f | sort)
  [ "${#family_gffs[@]}" -ge 2 ] || continue
  family="$(basename "$family_dir")"
  mkdir -p "$RESULTS/plasmid_panaroo/$family"
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$RESULTS/plasmid_panaroo/$family" \
    quay.io/biocontainers/panaroo:1.8.1--pyhdfd78af_0 \
    panaroo -i "${family_gffs[@]}" -o "$RESULTS/plasmid_panaroo/$family/panaroo_output" \
      --clean-mode sensitive --remove-invalid-genes --all_seq_in_graph \
      --len_dif_percent 0.95 --length_outlier_support_proportion 0.95 \
      --min_trailing_support 1 --trailing_recursive 1 \
      --min_edge_support_sv 0 --remove_by_consensus False \
      --edge_support_threshold 0 --threads 12
done
```

## 7. Call chromosome and plasmid variants with Snippy

Snippy uses the quality-trimmed paired short reads for both assembly
strategies; the long reads are used for hybrid assembly, not Snippy calling.

```bash
mkdir -p "$RESULTS/snippy/chromosome" "$RESULTS/snippy/plasmids"
chromosome_reference="$RESULTS/panaroo/chromosome/pan_genome_reference.fa"

for r1 in "$READS"/*_R1.fastq.gz; do
  sample="$(basename "$r1" _R1.fastq.gz)"
  sample_dir="$RESULTS/samples/$sample"
  docker run --rm --user "$DOCKER_USER" \
    -v "$RUN_ROOT:$RUN_ROOT" -w "$RESULTS/snippy/chromosome" \
    staphb/snippy:4.6.0 \
    snippy --cpus 12 \
      --outdir "$RESULTS/snippy/chromosome/${sample}_snippy_results" \
      --ref "$chromosome_reference" \
      --R1 "$sample_dir/fastp/${sample}_cleaned_R1.fastq.gz" \
      --R2 "$sample_dir/fastp/${sample}_cleaned_R2.fastq.gz" \
      --prefix "$sample" --mapqual 60 --basequal 13 \
      --mincov 10 --minfrac 0.9 --minqual 100 --maxsoft 100 \
      --report --unmapped \
      --fbopt "--min-alternate-count 2 --min-alternate-fraction 0.2 --min-base-quality 13" \
      --force
done

for family_dir in "$RESULTS"/plasmids/families/*; do
  [ -d "$family_dir" ] || continue
  family="$(basename "$family_dir")"
  reference="$RESULTS/plasmid_panaroo/$family/panaroo_output/pan_genome_reference.fa"
  [ -s "$reference" ] || continue
  for gff in "$family_dir"/*.gff3; do
    [ -f "$gff" ] || continue
    plasmid="$(basename "$gff" .gff3)"
    sample="${plasmid%%__*}"
    sample_dir="$RESULTS/samples/$sample"
    outdir="$RESULTS/snippy/plasmids/$family/${plasmid}_snippy_results"
    mkdir -p "$outdir"
    docker run --rm --user "$DOCKER_USER" \
      -v "$RUN_ROOT:$RUN_ROOT" -w "$RESULTS/snippy/plasmids/$family" \
      staphb/snippy:4.6.0 \
      snippy --cpus 12 \
        --outdir "$outdir" \
        --ref "$reference" \
        --R1 "$sample_dir/fastp/${sample}_cleaned_R1.fastq.gz" \
        --R2 "$sample_dir/fastp/${sample}_cleaned_R2.fastq.gz" \
        --prefix "$plasmid" --mapqual 60 --basequal 13 \
        --mincov 10 --minfrac 0.9 --minqual 100 --maxsoft 100 \
        --report --unmapped \
        --fbopt "--min-alternate-count 2 --min-alternate-fraction 0.2 --min-base-quality 13" \
        --force
  done
done
```

## 8. Inspect and clean temporary files

```bash
find "$RESULTS" -maxdepth 4 -type f -print
```

After confirming the outputs, remove only temporary process data you no longer
need. Keep the results, input reads, and databases:

```bash
rm -rf "$SRA_CACHE" "$SRA_TMP" "$RESULTS/mge/sequences"
```
