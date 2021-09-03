# This is the same code that vep uses in its wrapper, but we remove the plugins imput  due to 
# we are not using them.

from pathlib import Path
from snakemake.shell import shell


def get_only_child_dir(path):
    children = [child for child in path.iterdir() if child.is_dir()]
    assert (
        len(children) == 1
    ), "Invalid VEP cache directory, only a single entry is allowed, make sure that cache was created with the snakemake VEP cache wrapper"
    return children[0]


extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

fork = "--fork {}".format(snakemake.threads) if snakemake.threads > 1 else ""
stats = snakemake.output.stats
cache = snakemake.input.get("cache", "")


if snakemake.output.calls.endswith(".vcf.gz"):
    fmt = "z"
elif snakemake.output.calls.endswith(".bcf"):
    fmt = "b"
else:
    fmt = "v"

fasta = snakemake.input.get("fasta", "")
if fasta:
    fasta = "--fasta {}".format(fasta)

gff = snakemake.input.get("gff", "")
if gff:
    gff = "--gff {}".format(gff)

if cache:
    entrypath = get_only_child_dir(get_only_child_dir(Path(cache)))
    species = entrypath.parent.name
    release, build = entrypath.name.split("_")
    cache = (
        "--offline --cache --dir_cache {cache} --cache_version {release} --species {species} --assembly {build}"
    ).format(cache=cache, release=release, build=build, species=species)

shell(
    "(vep {extra} {fork} "
    "--format vcf "
    "{cache} "
    "{gff} "
    "{fasta} "
    "--output_file STDOUT "
    "--stats_file {stats} "
    "-i {snakemake.input.calls} "
    "-o {snakemake.output.calls}) {log}"
)
