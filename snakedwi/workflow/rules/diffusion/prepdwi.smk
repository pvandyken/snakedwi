
wildcard_constraints:
    shell="[0-9]+",


rule import_dwi:
    input:
        dwi=re.sub(".nii.gz", ".{ext}", input_path["dwi"]),
    output:
        dwi=temp(bids(
            root=work,
            suffix="dwi.{ext,nii.gz|bval|bvec|json}",
            datatype="dwi",
            **input_wildcards["dwi"]
        )),
    group:
        "subj"
    shell:
        "cp {input.dwi} {output.dwi}"


rule dwidenoise:
    input:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                datatype="dwi",
                **input_wildcards["dwi"],
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    output:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                desc="denoise",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="denoise.log", **input_wildcards["dwi"]),
    group:
        "subj"
    shell:
        "dwidenoise {input[0]} {output[0]} 2> {log} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"


def get_degibbs_inputs(wildcards):
    # Denoise as input if at least 30 dirs or not skipped
    import numpy as np

    in_dwi_bval = re.sub(".nii.gz", ".bval", input_path["dwi"].format(**wildcards))
    bvals = np.loadtxt(in_dwi_bval)

    if bvals.size < 30 or config["skip_denoise"]:
        prefix = bids(root=work, suffix="dwi", datatype="dwi", **wildcards)
    else:
        prefix = bids(
            root=work, suffix="dwi", datatype="dwi", desc="denoise", **wildcards
        )
    return multiext(prefix, ".nii.gz", ".bvec", ".bval", ".json")


rule mrdegibbs:
    input:
        get_degibbs_inputs,
    output:
        multiext(
            bids(
                root=work,
                suffix="dwi",
                datatype="dwi",
                desc="degibbs",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="degibbs.log", **input_wildcards["dwi"]),
    group:
        "subj"
    shell:
        "mrdegibbs {input[0]} {output[0]} 2> {log} && "
        "cp {input[1]} {output[1]} && "
        "cp {input[2]} {output[2]} && "
        "cp {input[3]} {output[3]}"


def _get_ref_run(wildcards):
    filtered = inputs["dwi"].filter(
        subject=wildcards.get("subject"), session=wildcards.get("session")
    )
    sizes = [
        os.path.getsize(p) for p in expand(filtered.path, zip, **filtered.zip_lists)
    ]
    # adapted from https://stackoverflow.com/a/3382369
    largest_ix = sorted(range(len(sizes)), key=sizes.__getitem__)[-1]

    ref = expand(
        bids(
            root=work,
            suffix="b0.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
        zip,
        **filtered.zip_lists
    )[largest_ix]
    return ref


rule rigid_align_runs:
    input:
        ref=_get_ref_run,
        b0=bids(
            root=work,
            suffix="b0.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **input_wildcards["dwi"]
        ),
        data=multiext(
            bids(
                root=work,
                suffix="dwi",
                desc="degibbs",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ),
    params:
        do_txf=lambda wildcards, input: (
            "true" if input["ref"] != input["b0"] else "false"
        )
    output:
        b0=temp(bids(
            root=work,
            suffix="b0reg.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **input_wildcards["dwi"]
        )),
        data=temp(multiext(
            bids(
                root=work,
                suffix="dwireg",
                desc="degibbs",
                datatype="dwi",
                **input_wildcards["dwi"]
            ),
            ".nii.gz",
            ".bvec",
            ".bval",
            ".json",
        ))
    shadow: 'minimal'
    container:
        config["singularity"]["mrtrix"]
    group:
        "subj"
    shell:
        """
        if {params.do_txf}; then
            mrregister -type rigid {input.b0} {input.ref} -rigid txf.mat
            mrtransform {input.b0} {output.b0} -template {input.ref} -linear txf.mat
            mrconvert {input.data[0]} -fslgrad {input.data[1]} {input.data[2]} data.mif
            mrtransform data.mif resliced.mif -template {input.ref} -linear txf.mat 
            mrconvert resliced.mif {output.data[0]} \\
                -export_grad_fsl {output.data[1]} {output.data[2]}
        else
            cp {input.b0} {output.b0}
            cp {input.data[0]} {output.data[0]}
            cp {input.data[1]} {output.data[1]}
            cp {input.data[2]} {output.data[2]}
        fi
        cp {input.data[3]} {output.data[3]}
        """


def get_concat_or_cp_cmd(wildcards, input, output):
    """Concatenate (if multiple inputs) or copy"""
    if len(input) > 1:
        cmd = f"mrcat"
    elif len(input) == 1:
        cmd = f"cp"
    else:
        # no inputs
        cmd = None
    return cmd


rule concat_degibbs_dwi:
    input:
        dwi_niis=lambda wildcards: get_dwi_indices(
            inputs["dwi"].filter(**wildcards).expand(
                bids(
                    root=work,
                    suffix="dwireg.nii.gz",
                    desc="degibbs",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
            ),
            wildcards,
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        dwi_concat=bids(
            root=work,
            suffix="dwireg.nii.gz",
            desc="degibbs",
            datatype="dwi",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="concat_degibbs_dwi.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} {input} {output} 2> {log}"


rule concat_runs_bvec:
    input:
        bv_files=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwireg.bvec",
                    desc="{{desc}}",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        out_fname=bids(
            root=work,
            suffix="dwi.bvec",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/concat_bv.py"


rule concat_runs_bval:
    input:
        bv_files=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwireg.bval",
                    desc="{{desc}}",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        out_fname=bids(
            root=work,
            suffix="dwi.bval",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/concat_bv.py"


# Combine multiple json from multiple scans (currently only copying first)
rule concat_runs_json:
    input:
        jsons=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="dwi.json",
                    desc="{{desc}}",
                    datatype="dwi",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        json=bids(
            root=work,
            suffix="dwi.json",
            desc="{desc}",
            datatype="dwi",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cp {input.jsons[0]} {output.json}"


rule get_shells_from_bvals:
    input:
        bval="{dwi_prefix}.bval",
    output:
        json="{dwi_prefix}.shells.json",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shells_from_bvals.py"


# Write 4D dwi_file with average shells
rule get_shell_avgs:
    input:
        dwi="{dwi_prefix}.nii.gz",
        shells="{dwi_prefix}.shells.json",
    output:
        avgshells="{dwi_prefix}.avgshells.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shell_avgs.py"


# Extract individual shell (e.g. B0) before rigid alignment (output used to calculate alignemnt)
rule get_shell_avg:
    input:
        dwi="{dwi_prefix}_dwi.nii.gz",
        shells="{dwi_prefix}_dwi.shells.json",
    params:
        bval="{shell}",
    output:
        avgshell="{dwi_prefix}_b{shell}.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shell_avg.py"


# Extract vols from particular shell (e.g. B0), after rigid alignment
rule get_shell_vols:
    input:
        dwi="{dwi_prefix}_dwireg.nii.gz",
        shells="{dwi_prefix}_dwi.shells.json",
    params:
        bval="{shell}",
    output:
        shell_vols="{dwi_prefix}_b{shell}s.nii.gz",
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_shell_vols.py"


# now have nii with just the b0's, want to create the topup phase-encoding text files for each one:
rule get_phase_encode_txt:
    input:
        bzero_nii=bids(
            root=work,
            suffix="b0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
        json=bids(
            root=work,
            suffix="dwi.json",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
    output:
        phenc_txt=bids(
            root=work,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **input_wildcards["dwi"]
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/diffusion/get_phase_encode_txt.py"


rule concat_phase_encode_txt:
    input:
        phenc_txts=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="phenc.txt",
                    datatype="dwi",
                    desc="degibbs",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    output:
        phenc_concat=bids(
            root=work,
            suffix="phenc.txt",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    group:
        "subj"
    shell:
        "cat {input} > {output}"


rule concat_bzeros:
    input:
        bzero_niis=lambda wildcards: get_dwi_indices(
            expand(
                bids(
                    root=work,
                    suffix="b0reg.nii.gz",
                    datatype="dwi",
                    desc="degibbs",
                    **input_wildcards["dwi"]
                ),
                zip,
                **filter_list(input_zip_lists["dwi"], wildcards)
            ),
            wildcards,
        ),
    params:
        cmd=get_concat_or_cp_cmd,
    output:
        bzero_concat=bids(
            root=work,
            suffix="concatb0.nii.gz",
            datatype="dwi",
            desc="degibbs",
            **subj_wildcards
        ),
    container:
        config["singularity"]["mrtrix"]
    log:
        bids(root="logs", suffix="concat_bzeros.log", **subj_wildcards),
    group:
        "subj"
    shell:
        "{params.cmd} {input} {output} 2> {log}"


def get_b0_mask():
    # Method options
    methods = {
        "b0_BET": "bet_from-b0",
        "b0_SyN": f"b0SyN_from-{config['template']}",
        "b0_synthstrip": "synthstrip_from-dwirefb0",
    }

    # Get BIDS name of file
    return bids(
        root=work,
        suffix="mask.nii.gz",
        desc="brain",
        method=methods.get(config["masking_method"]),
        datatype="dwi",
        **subj_wildcards,
    )


rule qc_b0_brainmask:
    input:
        img=rules.cp_dwi_ref.output.dwi_ref,
        seg=get_b0_mask(),
    output:
        png=report(
            bids(
                root=root,
                datatype="qc",
                suffix="mask.png",
                desc="brain",
                **subj_wildcards
            ),
            caption="../report/brainmask_dwi.rst",
            category="Brainmask",
        ),
        html=bids(
            root=root,
            datatype="qc",
            suffix="mask.html",
            desc="brain",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../../scripts/qc/vis_qc_dseg.py"


rule compile_qc_b0_brainmask_manifest:
    input:
        mask_qc=expand(rules.qc_b0_brainmask.output.png, zip, **subj_zip_list),
    output:
        os.path.join(qc, "data", "mask.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "title": "DWI Brainmask",
                    "images": sorted(input["mask_qc"]),
                },
                f,
            )
