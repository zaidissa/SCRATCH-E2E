/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    asBool — read a flag that may have come from the command line
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    A param declared `skip_x = false` in nextflow.config is a Boolean, but the
    same param supplied as `--skip_x false` arrives as the STRING "false", and
    Nextflow does not coerce it to the declared type. In Groovy every non-empty
    string is truthy, so

        if (!params.skip_x)          // "false" -> !"false" -> false

    takes the same branch for "true" and for "false": passing the flag at all
    disabled the stage, whatever value you gave it. Verified against this
    pipeline's own config:

        default            skip_celltypist = true  (Boolean)
        --skip_celltypist false             = false (String)  -> !v == false

    which is why --skip_celltypist false could never turn CellTypist back on.

    Every flag read in a conditional goes through here. Config files cannot
    declare functions (Nextflow 26.x rejects it), so `conf/*.config` inlines the
    same expression instead.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def asBool(value) {
    if (value == null)            return false
    if (value instanceof Boolean) return value
    return value.toString().trim().toLowerCase() in ['true', 'yes', 'y', 'on', '1']
}
