
rule qc_v1_rgb:
    input:
        fa=rules.dtifit_resampled_t1w.output.out_fa,
        v1=rules.dtifit_resampled_t1w.output.out_v1,
    output:
        temp(bids(
            root=root,
            datatype="qc",
            desc="preproc",
            suffix="v1rgb.nii.gz",
            **subj_wildcards
        )),
    container:
        config["singularity"]["mrtrix"]
    shell:
        """
        mrcalc {input.fa} {input.v1} -multiply -abs {output}
        """

rule qc_fa:
    input:
        img=rules.qc_v1_rgb.output[0]
    output:
        png=bids(
            root=root,
            datatype="qc",
            suffix="fa.png",
            desc="preproc",
            **subj_wildcards
        ),
    group:
        "subj"
    container:
        config["singularity"]["python"]
    script:
        "../scripts/qc/vis_qc_rgb.py"


rule compile_qc_fa:
    input:
        imgs=expand(rules.qc_fa.output.png, zip, **subj_zip_list),
    output:
        os.path.join(qc, "data", "fa.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "title": "FA",
                    "images": sorted(input["imgs"]),
                },
                f,
            )

rule qc:
    input:
        mask_qc=rules.compile_qc_b0_brainmask_manifest.output[0],
        reg_qc=rules.compile_qc_reg_dwi_t1_manifest.output[0],
        fa_qc=rules.compile_qc_fa.output[0],
    output:
        os.path.join(qc, "data.json"),
    run:
        with open(output[0], "w") as f:
            json.dump(
                {
                    "mask": json.loads(Path(input["mask_qc"]).read_text()),
                    "reg": json.loads(Path(input["reg_qc"]).read_text()),
                    "fa": json.loads(Path(input["fa_qc"]).read_text()),
                },
                f,
            )


_qc_app = os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz")


def _get_tar_contents(file):
    try:
        return [
            p
            for p in sp.check_output(["tar", "-tf", _qc_app]).decode().splitlines()
            if p[-1] != "/"
        ]
    except sp.CalledProcessError as err:
        raise Exception("Unable to find qc-app.tar.gz...") from err


rule unpack_qc_app:
    input:
        os.path.join(workflow.basedir, "..", "resources", "qc-app.tar.gz"),
    output:
        _get_tar_contents(_qc_app),
    shell:
        "tar -xvzf {input}"
