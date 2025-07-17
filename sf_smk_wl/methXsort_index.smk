
graft_fa = config["graft_fa"]
host_fa = config["host_fa"]
outdir = os.path.abspath(config["outdir"])
xengsort_index_prefix = config.get("xengsort_index_prefix", None)

methxsort_path = config["methxsort_path"]

rule all:
    input:
        os.path.join(outdir, "ref_idx/xengsort_index.log"),

rule convert_ref:
    input:
        graft_fa = graft_fa,
        host_fa = host_fa,
    output:
        graft_fa_cvt = os.path.join(outdir, "ref_idx/graft_cvt.fa"),
        host_fa_cvt = os.path.join(outdir, "ref_idx/host_cvt.fa"),
    shell:
        """
python {methxsort_path} convert-ref \
    {input.graft_fa} --out {output.graft_fa_cvt} 
python {methxsort_path} convert-ref \
    {input.host_fa} --out {output.host_fa_cvt}
"""

rule xengsort_idx:
    input:
        graft_fa = rules.convert_ref.output.graft_fa_cvt,
        host_fa = rules.convert_ref.output.host_fa_cvt,
    output:
        log = os.path.join(outdir, "ref_idx/xengsort_index.log"),
    params:
        xengsort_index_prefix = os.path.join(outdir, "ref_idx/xengsort_index"),
        build = 1

    shell:
        """
python {methxsort_path} xengsort-index \
    --host {input.host_fa} --graft {input.graft_fa} \
    --index {params.xengsort_index_prefix} \
    -n 7_000_000_000 --fill 0.88 --statistics summary -k 25 \
    > {output.log} 2>&1
"""