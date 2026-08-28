/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CELLBENDER_TO_CELLRANGER — re-emit CellBender counts in cellranger h5 layout

    Seurat's Read10X_h5 treats every ROOT-LEVEL group in the file as a genome and
    reads "<genome>/data" from each one. A cellranger file has exactly one root
    group, `matrix`. CellBender adds three more — droplet_latents, global_latents
    and metadata — so the reader walks into droplet_latents, finds no `data`, and
    the QC render dies with

        An object with name droplet_latents/data does not exist in this group

    The `matrix` group itself is already the cellranger layout, so copying just
    that group yields a file every downstream reader accepts.

    This is a SEPARATE process rather than a few lines appended to CELLBENDER
    deliberately. Adding either a step or an output declaration to CELLBENDER
    invalidates its cache entry, and on this dataset that means re-running about
    thirteen hours of CPU inference to fix a thirty-second file rewrite.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process CELLBENDER_TO_CELLRANGER {

    tag "Reformatting ${sample_id} for cellranger readers"
    label 'process_single'

    input:
        tuple val(sample_id), path(csv_metrics), path(cb_h5)

    output:
        tuple val(sample_id), path(csv_metrics),
              path("cellbender_${sample_id}_cellranger_filtered.h5"), emit: corrected

    when:
        task.ext.when == null || task.ext.when

    script:
        """
        python3 <<'PYEOF'
import h5py

src = "${cb_h5}"
dst = "cellbender_${sample_id}_cellranger_filtered.h5"

with h5py.File(src, "r") as fi:
    if "matrix" not in fi:
        raise SystemExit(
            f"{src} has no 'matrix' group (root groups: {list(fi.keys())}). "
            "CellBender 0.3.x writes one; a different version may not.")
    with h5py.File(dst, "w") as fo:
        fi.copy("matrix", fo, name="matrix")
        for k, v in fi.attrs.items():      # chemistry_description, filetype, version
            fo.attrs[k] = v

# Assert the property that actually matters to the reader, rather than trusting
# the copy: exactly one root group, named matrix.
with h5py.File(dst, "r") as f:
    roots = list(f.keys())
    if roots != ["matrix"]:
        raise SystemExit(f"expected a single root group 'matrix', found {roots}")
    genes, cells = f["matrix"]["shape"][:]
    print(f"${sample_id}: {genes} genes x {cells} cells, cellranger layout")
PYEOF
        """

    stub:
        """
        touch cellbender_${sample_id}_cellranger_filtered.h5
        """
}
