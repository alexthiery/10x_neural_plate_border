static String build_debug_param_summary(params, monochrome_logs) {
    Set paramsKeySet = params.keySet()
    def summary = [:]
    paramsKeySet.each {
        if(it != "modules" && it != "genomes") {
            summary[it] = params.get(it)
        }
    }
    def output = summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    output += "\n" + dashed_line(monochrome_logs)
    return output
}

static String build_debug_scripts_summary(params, monochrome_logs) {
    Set paramsKeySet = params.keySet()
    def summary = [:]
    paramsKeySet.each {
            summary[it] = params.get(it)
    }
    def output = summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
    output += "\n" + dashed_line(monochrome_logs)
    return output
}

private static Map log_colours(Boolean monochrome_logs) {
    Map colorcodes = [:]
    colorcodes['reset']       = monochrome_logs ? '' : "\033[0m"
    colorcodes['dim']         = monochrome_logs ? '' : "\033[2m"
    colorcodes['black']       = monochrome_logs ? '' : "\033[0;30m"
    colorcodes['green']       = monochrome_logs ? '' : "\033[0;32m"
    colorcodes['yellow']      = monochrome_logs ? '' :  "\033[0;33m"
    colorcodes['yellow_bold'] = monochrome_logs ? '' : "\033[1;93m"
    colorcodes['blue']        = monochrome_logs ? '' : "\033[0;34m"
    colorcodes['purple']      = monochrome_logs ? '' : "\033[0;35m"
    colorcodes['cyan']        = monochrome_logs ? '' : "\033[0;36m"
    colorcodes['white']       = monochrome_logs ? '' : "\033[0;37m"
    colorcodes['red']         = monochrome_logs ? '' : "\033[1;91m"
    return colorcodes
}
static String dashed_line(monochrome_logs) {
    Map colors = log_colours(monochrome_logs)
    return "-${colors.dim}----------------------------------------------------${colors.reset}-"
}
