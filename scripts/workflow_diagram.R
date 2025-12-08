DiagrammeR::grViz(
  "digraph {

graph [layout = dot, rankdir = TB]

# define the global styles of the nodes. We can override these in box if we wish
node [shape = rectangle, style = filled, fillcolor = linen]

main [label = '00-main.R', shape = tab]

run_for [label = 'run_for', shape = diamond, fillcolor = white]

subgraph cluster_ab {
  label = 'Alberta'
  bgcolor = azure

  download_ab [label= '01-download-data.R', shape = tab]
  extract [label = '01a-extract-mdb.R', shape = tab]
  import [label= '01b-import-mdb-csv.R', shape = tab]
  dataprep_ab [label= '02a-Alberta-data-prep.R', shape = tab]
  explore_ab [label= '02b-Alberta-explore.R', shape = tab]
  analyses_ab [label= '02c-Alberta-analyses.R', shape = tab]
}

subgraph cluster_np {
  label = 'National Parks'
  bgcolor = lavender

  download_np [label= '01-download-data.R', shape = tab]
  analyses_np [label= '03-Jasper-analyses.R', shape = tab]
}

# edge definitions with the node IDs
main -> run_for -> {download_ab download_np}
download_ab -> extract -> import [style = dashed]
download_ab -> import -> dataprep_ab -> explore_ab -> analyses_ab
download_np -> analyses_np

}"
) |>
  DiagrammeRsvg::export_svg() |>
  charToRaw() |>
  rsvg::rsvg_png("workflow.png")
